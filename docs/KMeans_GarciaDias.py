#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os
import multiprocessing
import sys

def assign_array(clusters_tmp, fullset, i_worker, output):
	"""Assign vectors to classes minimizing euclidean distances.
	This line is run in parallel.
	
	Parameters
	----------
	clusters_tmp : 2D array-like, shape=(n_clusters, n_features)
		Array with the partial cluster centers (n_clusters, n_features).
	
	fullset : 2D array, shape=(n_samples, n_features)
		Array with all objects in to be classified.
	
	i_worker : int
		Identifier to the ith job. As it runs in parallel is important to index each result.

	output : object type
		Object keeping the results to each job.
	
        Returns
        -------
		The return in this function is implicit as output attributes.
		The output is a list (n_samples + 1) witch the first value correspond to ith job and other values are
		the classification of the objects in the sample.
	"""
	output.put([i_worker] + [(((clusters_tmp - line)**2).sum(axis=1)).argmin() for line in fullset])

def three_sigma_cut(n_clusters, assign, fullset_tmp, clusters, sclusters):
	"""3-sigma cut of the outliers of the calculation of the final cluster centres and mean errors.
	We just discard outliers in classes with more than 10 objects (npos >= 10), 10 been arbitrarily chosen.

	Parameters
	-----------
	n_clusters : int
		Number of clusters
	
	assign : 1D array, shape=(n_samples)
		Vector with the classification. For each spectrum a int number
		is assigned, relative to the class they belong.
	fullset_tmp : 2D array, shape=(n_samples, n_features)
		The data to pick seeds for. To avoid memory copy, the input data
		should be double precision (dtype=np.float64).
	clusters : array shape=(n_samples, n_features)
		Cluster centers. It is the avarege of the spectra in each class.
	sclusters : array shape=(n_samples, n_features)
		The standard deviation of each cluster.	
	Returns
	-------
	
	new_clusters : 2D array, shape=(n_samples, n_features)
		Gives each cluster center without tanking outliers in to account. It is the avarege of the spectra in each class.

	new_sclusters : 2D array, shape=(n_samples, n_features)
		The standard deviation of each cluster without tanking outliers in to account.	
	"""
        n_variables = clusters.shape[1]
	POS = np.array([], dtype='int')
        for i in xrange(n_clusters):
                pos = np.where(assign == i)[0]
                npos = len(pos)
                if (npos >= 10): #*Observation in reader refers to this line.
                        sigma = np.nanmean(sclusters[i,:])
                        res = np.ones(npos)
                        vec_temp = fullset_tmp[pos]
                        res = np.nanmean(abs(vec_temp - clusters[i]), axis=1)
                        ppos = np.where(res <= 3.0*sigma)[0]
			sys.stdout.flush()
                        print('%d outliers removed from class %d'% (len(pos) - len(ppos),i), end='\r')
                        pos = pos[ppos]
			POS = np.append(POS, pos)
		else:
			POS = np.append(POS, pos)
	new_clusters, new_sclusters = centroid_calc(n_clusters, assign[POS], fullset_tmp[POS])
        return new_clusters, new_sclusters

def centroid_calc(n_clusters, assign, fullset_tmp):
	"""Calculate centroides of the classes and the internal sigma in each dimension.

	Parameters
	-----------
	n_clusters: int
		Number of clusters
	
	assign: 1D array, shape=(n_samples)
		Vector with the classification. For each spectrum a int number
		is assigned, relative to the class they belong.

	fullset_tmp : 2D array, shape=(n_samples, n_features)
		Array with spectral data to be classified. Spectra are rows, columns are flux.

	Returns
	-------
	
	clusters : 2D array, shape=(n_samples, n_features)
		Gives each cluster center. It is the avarege of the spectra in
		each class.

	sclusters : 2D array, shape=(n_samples, n_features)
		The standard deviation of each cluster.	
	"""
	n_variables = fullset_tmp.shape[1]
        clusters = np.zeros((n_clusters, n_variables))
        sclusters = np.copy(clusters)
        zero_var = np.zeros(n_variables)
	EMPTY_CLASSES = ''
	N_empty = 0
        for i_class in xrange(n_clusters):
		clean = [assign == i_class]
                vec_temp = fullset_tmp[clean]
		if (len(vec_temp) != 0):
	                clusters[i_class]  = np.nanmean(vec_temp, axis=0)
			sclusters[i_class] = np.nanstd(vec_temp, axis=0)
		else:
			CLASS_STR = '%d, '%i_class
			EMPTY_CLASSES = ''.join((EMPTY_CLASSES, CLASS_STR))
			N_empty += 1
	if (N_empty == 1): print ('\nClass %s is empty!!!'%EMPTY_CLASSES[:-2])
	if (N_empty >  1): print ('\nClasses %s are empty!!!'%EMPTY_CLASSES[:-2])
        return clusters, sclusters

def K_MEANS_FOR_STARS(fullset, saveset='k_means_results', mask=[], n_clusters=20, centers0=0, RandState=None, N_CPU=1):
	"""Perform k-means classification of the APOGEE spectra.
	
	Parameters
	----------
	fullset : 2D array, shape=(n_samples, n_features)
		Array with spectral data to be classified. Spectra are rows, columns are flux.
	
	saveset : string, optional, default: 'k_means_results'
		Filename for saving outputs in npz-like file.
	
	mask : list, optinal, shape=(<n_features), default: []
		Empty list means all columns are used.
		List of column positions you want to use in the classification.
		Cuting the original vector can be very time cosuming.
		If you have a large dataset and want to run the code many times,
		you shoud reshape the fullset input instead of using a mask.
	
	n_clusters : int, optinal, default: 20
		Number of output clusters.
	
	centers0 : 2D array, shape=(n_clusters,n_features), default: 0
		Initial clusters. If centers0 == 0, then perform randomized initialization.

	RandState : float, optional, default: None
		Random state to be used in random initialization.

	N_CPU : int, optional, default: 1
		Number of jobs executed in parallel.

	Returns
	-------
	
	count: 1D array, shape=(n_clusters)
		Number of spectra in each cluster
	
	assign: 1D array, shape=(n_samples)
		Vector with the classification. For each spectrum a int number
		is assigned, relative to the class they belong.

	clusters: 2D array, shape=(n_clusters, n_samples)
		Gives each cluster center. It is the avarege of the spectra in
		each class.

	sclusters: 2D array, shape=(n_clusters, n_samples)
		The standard deviation of each cluster.
	"""
	soutput = saveset + '.npz'
	setc = fullset
	print('Starting kMeans in %d stars with %d classes. Executing %d parallel jobs. ' % (len(setc), n_clusters, N_CPU))
	vec_example = setc[0]
	n_variables = len(vec_example)
	n_samples = len(setc)
	n_iter = 40 #Maximum number of iterations
	#Executes k_means
	assign, clusters0 = K_MEANS7(fullset=setc, n_clusters=n_clusters, n_iterations=n_iter, mask=mask, centers0=centers0, RandState=RandState, N_CPU=N_CPU)
	#Start to look on the results
	clusters, sclusters = centroid_calc(n_clusters, assign, setc)
	#Discard outliers in the means account
	clusters, sclusters = three_sigma_cut(n_clusters, assign, setc, clusters, sclusters)
	count, CLASS_SORT, assign = map_classes(assign, len(clusters))
	print('Saving database in %s' % soutput)
	np.savez(soutput, assign=assign, clusters=clusters,
			sclusters=sclusters, n_clusters=n_clusters, count=count, centers0=centers0)
	return clusters, assign

def K_MEANS7(fullset, n_clusters=10, n_iterations=2, mask=0, centers0=0, RandState=None, N_CPU=1):
	"""Cluster centers determination
	
	Parameters
	----------
	fullset: array or sparse matrix, shape (n_samples, n_features)
		spectra to be classified

	n_clusters: int, optional, default = 20
		Number of cluster centers.

	n_iterations: int, optional, default = 10
		Use this keyword to specify the maximum
                number of iterations in computing the cluster
                centers.

	mask : list, shape=(<n_features), optinal, default: []
		Empty list means all columns are used.
		List of column positions you want to use in the classification.
		Cuting the original vector can be very time cosuming.
		If you have a large dataset and want to run the code many times,
		you shoud reshape the fullset input instead of using a mask.
	
	n_clusters : int, optinal, default: 20
		Number of output clusters.
	
	centers0 : 2D array, shape=(n_clusters,n_features), default: 0
		Initial clusters. If centers0 == 0, then perform randomized initialization.

	RandState : float, optional, default: None
		Random state to be used in random initialization.

	N_CPU : int, optional, default: 1
		Number of jobs executed in parallel.

	Returns
	-------
	clusters : 2D array, shape=(n_clusters, n_features)
		Gives the cluster centers.
		The clusters are sorted according to the
                number of members, 0 being the most numerous,
                1 the second most numerous, and so on.

	assign_pos : vector, shape=(n_clusters, n_features)
		Pointer where each element of array is
		assigned to one of the cluster_centers.

	Notes
	-----
		This functions just states what are the clusters centers used in
		the final classification. The function applies the k-mean classification
		with a limited number of centers and iteration so we can quickly define
		what are the moust robust centers to apply the classification.
	"""

	vec_example = fullset[0]
	if (len(mask) != 0): vec_example = vec_example[mask]
	dimension = [len(vec_example), len(fullset)]
	if (len(dimension) != 2):
		print('Input array must be two-dimensional. \
			It has %d\t%d' % (dimension[0], dimension[1])  )
	n_variables = dimension[0]
	n_samples = dimension[1]
	ZERO_VAR_CLUSTER = np.zeros((n_clusters, n_variables), dtype='Float32')
	if (len(np.shape(centers0)) == 0):
		print('Starting my own initialization')

		www0 = np.copy(ZERO_VAR_CLUSTER)
		n_iter = 2 #max number of iterations
		n_cls = 10 #number of classes. Keep in mind that you just retain 1.
		ini_cls = 0
		npos = n_clusters + 1
		
		while (npos > 0 and ini_cls < n_clusters):
			if (ini_cls == 0):
				count, assign_pos, clusters_tmp, n_iter = ACTUAL_K_MEANS(fullset,
										n_clusters=n_cls, 	
										n_iterations=n_iter,
										mask = mask,
										RandState=RandState,
										N_CPU=N_CPU)
			else:
				count, assign_pos, clusters_tmp, n_iter = ACTUAL_K_MEANS(subarr,
										n_clusters=n_cls, 	
										n_iterations=n_iter,
										mask = mask,
										RandState=RandState,
										N_CPU=N_CPU)
			www0[ini_cls,:] = clusters_tmp[0,:]
			pos = np.where(assign_pos != 0)[0]
			npos = len(pos)
			if (npos != 0):
				if (ini_cls == 0):
					subarr = fullset
				else:
					subarr = fullset[pos]
			ini_cls = ini_cls + 1
			if (npos < 0.05*len(fullset)):
				for j in xrange(ini_cls, n_clusters, 1):
					try: vec_temp = subarr[np.random.randint(len(subarr))]
					except UnboundLocalError: vec_temp = fullset[np.random.randint(len(fullset))]
					#vec_temp = np.load( name)['norm']
					if (len(mask) != 0): vec_temp = vec_temp[mask]
					www0[j,:] = vec_temp
				ini_cls = n_clusters
		clusters = np.asarray(www0)
		n_clusters = np.shape(clusters)[0]
	else:
		clusters = np.zeros((n_variables, n_clusters), dtype='Float32')
		for i_name, name in enumerate(centers0):
			vec_temp = np.load( name)['norm']
			if (len(mask) != 0): vec_temp = vec_temp[mask]
			clusters[:,i_name] = vec_temp
	clusters0 = np.copy(clusters)
	count, assign_pos, clusters_tmp, n_iterations = ACTUAL_K_MEANS(fullset,
								clusters0=clusters,
								n_iterations=n_iterations,
								n_clusters=n_clusters,
								mask=mask,
								RandState=RandState,
								N_CPU=N_CPU)
	
	return assign_pos, clusters0

def ACTUAL_K_MEANS(fullset, clusters0=0, n_clusters=20, n_iterations=10, mask=0, RandState=None, N_CPU=1):
	"""Actual classification

	Note
	----
		This contains the actual classification procedure,
                which must be separated from the actual k_means
                since it has to be invoked many times for
                initialization.

		Inputs and outputs defined in K_means.

                Random initialization in case no input cluster center is given.

	Parameters
	----------
	fullset : 2D arrayr, shaper=(n_samples, n_features)
		spectra to be classified

	clusters0 : 2D array, shape=(n_clusters, n_features), optional, default: 0
		Gives the initial cluster centers, if 0 is assigned then
		centers are randomly assigned.

	n_clusters : int, optional, default: 20
		Number of cluster centers.

	n_iterations : int, optional, default: 10
		Use this keyword to specify the maximum
                number of iterations in computing the cluster
                centers.

	mask : list, shape=(<n_features), optinal, default: []
		Empty list means all columns are used.
		List of column positions you want to use in the classification.
		Cuting the original vector can be very time cosuming.
		If you have a large dataset and want to run the code many times,
		you shoud reshape the fullset input instead of using a mask.
	
	RandState : float, optional, default: None
		Random state to be used in random initialization.

	N_CPU : int, optional, default: 1
		Number of jobs executed in parallel.
	
	Returns
	-------
	count : 1D array, shape=(n_clusters)
		The number of members in each cluster.

	assign_pos :  1D array, shape=(n_samples)
		Pointer where each element of array is
		assigned to one of the cluster_centers.
	
	clusters_tmp : 2D array, shape=(n_clusters, n_samples)
		Gives the cluster centers.
		The clusters are sorted according to the
                number of members, 0 being the most numerous,
                1 the second most numerous, and so on.
	"""

	vec_example = fullset[0]
	if (len(mask) != 0): vec_example = vec_example[mask]
	
	dimension = [len(vec_example), len(fullset)]
	if (len(dimension) != 2 ): print('Input array must be two-dimensional')
	n_variables = dimension[0]
	n_samples = dimension[1]

	workcol = np.ones(n_clusters, dtype='Float32')
	workrow = np.ones(n_variables, dtype='Float32')
	assign_pos = np.zeros(n_samples, dtype='int')
	count = np.zeros(n_clusters, dtype='Float32')
	
	ZERO_VAR_CLUSTER = np.zeros((n_clusters, n_variables))
        if (len(np.shape(clusters0)) == 0):
		clusters_tmp_spec = []
		np.random.seed(RandState)
		for rand_value in np.random.randint(n_samples, size=n_clusters):
			if (len(mask) != 0):
				clusters_tmp_spec.append(fullset[rand_value][mask])
				np.random.seed(RandState)
			else:
				clusters_tmp_spec.append(fullset[rand_value])
				np.random.seed(RandState)
		clusters_tmp = np.asarray(clusters_tmp_spec)
        else:
        	print('\nImposing starting cluster centers!')
        	ssize = np.shape(clusters0)
        	if (ssize[1] != n_variables or ssize[0] != n_clusters):
        		print('Fix up this problem - centers have the wrong shape')
			quit()
        	else: clusters_tmp = np.copy(clusters0)

	#Iteration loop

	k_iter = 0
	caca = 1
	sample_list = range(n_samples)
	cluster_list = range(n_clusters)
	array_split = np.array_split(fullset, N_CPU)
	while (caca > 0.01 and k_iter < n_iterations):
		k_iter += 1
		assign_pos0 = np.copy(assign_pos)
		assign_pos = []
		clusters_tmp_tmp =  np.copy(ZERO_VAR_CLUSTER)
		output = multiprocessing.Queue()
		processes = [multiprocessing.Process(target=assign_array,
				args=(clusters_tmp, array_split[i_worker], i_worker, output)) for i_worker in range(N_CPU)]
		assign_tmp = []
		for p in processes: p.start()
		for p in processes: assign_tmp.append(output.get())
		for p in processes: p.join()
		assign_order = []
		for work in assign_tmp:
			assign_order.append(work[0])
		order = [i[0] for i in sorted(enumerate(assign_order), key=lambda x:x[1])]
		for i_order in order: assign_pos += assign_tmp[i_order][1:]

		assign_pos = np.asarray(assign_pos)
		for i_sample, line in enumerate(fullset):
			clusters_tmp_tmp[assign_pos[i_sample],:] += line
                for i in cluster_list:
                	pos = np.where(assign_pos == i)[0]
                	npos = len(pos)
                	count[i] = npos
                	if (npos > 0):
                		clusters_tmp[i,:] = clusters_tmp_tmp[i,:]/(1.*npos)
                	else: clusters_tmp[i,:] = clusters_tmp_tmp[i,:]
		pepe = np.where(assign_pos[:] - assign_pos0[:] != 0)[0]
		ncaca = len(pepe)
		caca = ncaca/(1.*n_samples)
		if (k_iter > 1): sys.stdout.flush()
		print("Iteration %2d, %5.1f%% of the classifications have been changed." % (k_iter, caca*100.), end='\r')

	n_iterations = k_iter
	count, CLASS_SORT, assign_pos = map_classes(assign_pos, n_clusters)
	clusters_tmp = clusters_tmp[CLASS_SORT,:]
	return count, assign_pos, clusters_tmp, n_iterations

def map_classes(assign,n_clusters):
       	'''Sort classes by number of objects
	
	Parameters
	----------
	
	assign : 1D array, shape=(n_samples)
		vector of assigned clusters for each spectra.
		
	n_clusters : int
		Number of clusters.
	
	Returns
	-------
	
	Count_CLASS : 1D array, shape=(n_clusters)
		Sorted number of members.
	
	CLASS_SORT : 1D array, shape=(n_clusters)
		Cluster indexes sorted by number of members.

	new_assign : 1D array, shape=(n_samples)
		New ordered assign of assigned clusters for each spectra.
	'''

 	Count_CLASS = np.zeros(n_clusters)
        new_assign = np.ones(len(assign))*9999
        POS = []
        for i in xrange(n_clusters):
                pos = np.where(assign == i)[0]
                Count_CLASS[i] = len(pos)
                POS.append(pos)
        CLASS_SORT = [order[0] for order in sorted(enumerate(Count_CLASS),key=lambda i:i[1])][::-1]
        POS = np.array(POS)
        for i_class, filt in enumerate(POS[CLASS_SORT]):
                new_assign[filt] = np.ones(len(filt))*i_class
        return Count_CLASS[CLASS_SORT], CLASS_SORT, new_assign

if __name__ == '__main__':
	print ('This code is running in pc:', end='\n')
	HOSTPC = os.system('hostname')
	NC = 50
	i = 42
	fullset = np.load('/net/tarea/scratch/Rafael/phd/apogee/python/FULL_masked.npy')[:1000]
	SET_RESULTS0 = K_MEANS_FOR_STARS(fullset,
					mask = [],
					saveset='/net/tarea/scratch/Rafael/phd/final_code/1000_sect_NC_'+str(NC)+'_Run_' + str(i),
					n_clusters = NC,
					RandState=i,
					N_CPU=6)
