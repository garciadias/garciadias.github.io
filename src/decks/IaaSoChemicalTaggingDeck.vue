<script setup>
import DeckBrand from '@/components/DeckBrand.vue';
import RevealDeck from '@/components/RevealDeck.vue';

// IAA-SO School on AI/ML in Astronomy 2026 · Unsupervised Learning pillar.
// A ~90-minute lecture. Narrative: stars born together share a chemical
// fingerprint, and clusters dissolve over time — so Galactic archaeology
// wants to *reconstruct* them from abundances alone (chemical tagging,
// Freeman & Bland-Hawthorn 2002; Kos et al. 2017). That is a clustering
// problem in a high-D, noisy space, so the talk is a tour of the clustering
// toolbox, each algorithm chosen to fix the failure mode of the last:
//
//   K-means (centroid, needs K) → KNN (the nearest-neighbour primitive) →
//   DBSCAN (density,
//   arbitrary shape) → HDBSCAN* (density at all scales) → PLSCAN
//   (persistence, no min-cluster-size) → t-SNE (project then see) → UMAP
//   (graph embedding, global) → EVoC (fuses UMAP + HDBSCAN* + PLSCAN).
//
// Source material: the Garcia-Dias et al. 2020 "Clustering analysis" chapter
// (K-means, KNN, DBSCAN, SSE/silhouette/dip/homogeneity), Bot et al. 2025
// (PLSCAN, persistence-based multiscale density clustering), and Kos et al.
// 2017 (t-SNE chemical tagging) — then the benchmark we run in this repo:
// t-SNE vs UMAP vs EVoC on APOGEE DR17 + Gaia EDR3, scored against kinematic
// ground truth.
//
// Style follows the PyData lightning deck (same RevealDeck + DeckBrand
// skeleton, same panel/checklist/crosslist vocabulary); the argument arc
// echoes the BigDataLDN talk ("why do some ML models fail" → no free lunch
// → know your data → never quit thinking), now pointed at stars.
const asset = (name) => `${import.meta.env.BASE_URL}presentations/iaa-so-chemical-tagging-2026/${name}`
</script>

<template>
  <RevealDeck :options="{ center: true }">
    <template #chrome>
      <DeckBrand :qr="asset('presentation.png')" :show-logos="false" />
    </template>

    <!-- 0 · Title -->
    <section class="title-slide center">
      <div class="eyebrow">IAA-SO School on AI/ML in Astronomy 2026 · Unsupervised Learning pillar</div>
      <h1>Chemical tagging: finding lost star clusters</h1>
      <p class="subtitle">
        Stars born together share a chemical fingerprint. Clusters dissolve, chemistry
        doesn't — so we reconstruct them from abundances alone, working our way from
        <strong>K-means to EVoC</strong>.
      </p>
      <p class="venue-note">Rafael Garcia-Dias · IAA-SO 2026 · ~90-minute lecture + hands-on</p>
      <p class="small muted">
        K-means → KNN → DBSCAN → HDBSCAN* → PLSCAN → t-SNE → UMAP → EVoC
      </p>
      <aside class="notes">
        (~2 min) One-sentence hook: a cluster's chemistry is a fossil that survives long after
        the cluster itself has scattered into the field. Frame the 90 minutes up front: we tour
        the clustering toolbox one algorithm at a time — each chosen to fix the last one's
        failure — and land on EVoC, the method that does embedding + clustering + scale
        selection in one pass. Then we benchmark it on real stars.
      </aside>
    </section>

    <!-- 1 · The science question -->
    <section>
      <div class="eyebrow">The science question</div>
      <h2>Stars are born together, then drift apart</h2>
      <div class="cols" style="--n: 3; margin-top: 0.3em">
        <div class="panel">
          <h3>The premise</h3>
          <p class="small">
            A cluster forms from one well-mixed cloud, so its members are chemically
            homogeneous to better than 0.1 dex — an abundance pattern is a
            <strong>fingerprint of common origin</strong>.
          </p>
        </div>
        <div class="panel">
          <h3>The dream</h3>
          <p class="small">
            <em>Galactic archaeology</em>: use that fingerprint to reconstruct clusters that
            dissolved billions of years ago, and to recover members flung far from the cluster
            centre.
          </p>
        </div>
        <div class="panel flip">
          <h3>Why it is a clustering problem</h3>
          <p class="small">
            Nothing tells us how many clusters there are, or which star belongs to which. We
            only have a cloud of points that should hide
            <strong>unlabelled overdensities</strong>.
          </p>
        </div>
      </div>
      <p class="small muted center" style="margin-top: 0.5em">
        <a href="https://arxiv.org/abs/1709.00794" target="_blank" rel="noopener">Kos et al. 2017</a>,
        after <a href="https://ui.adsabs.harvard.edu/abs/2002ARA%26A..40..487F/abstract" target="_blank" rel="noopener">Freeman &amp; Bland-Hawthorn 2002</a>.
      </p>
      <aside class="notes">
        (~2 min) Motivation first, no ML yet. Say the definition out loud: chemical tagging is
        using the abundances measured in stellar atmospheres to reconstruct chemically
        homogeneous clusters that have already dispersed. Then land the third panel, because it
        is the hinge of the whole lecture: we have no labels, no cluster count and no shapes —
        only points. That is the definition of a clustering problem, and everything from here on
        is about which clustering algorithm survives contact with this particular cloud.
      </aside>
    </section>

    <!-- 2 · The data -->
    <section>
      <div class="eyebrow">The data</div>
      <h2>A 16-dimensional chemical space (C-space)</h2>
      <div class="fig-split" style="--cols: 1fr 1fr; margin-top: 0.3em; align-items: start">
        <div>
          <p class="small">
            One star → one point in <strong>C-space</strong>, the space of its elemental
            abundances relative to iron. <strong>APOGEE DR17</strong> supplies the chemistry;
            <strong>Gaia EDR3</strong> astrometry is held back as ground truth.
          </p>
          <ul class="dotlist small">
            <li>16 abundance dimensions per star</li>
            <li>~183 000 stars after quality cuts (SNR ≥ 100)</li>
            <li>Each element standardised before any distance</li>
          </ul>
          <p class="small muted" style="margin-top: 0.4em">
            C, N, O, Na, Mg, Al, Si, S, K, Ca, Ti, V, Cr, Mn and Ni as [X/Fe], plus [Fe/H]
            itself. Rescaling every element to zero median and unit σ stops the one with the
            widest scatter from deciding every distance.
          </p>
        </div>
        <div class="figure">
          <img
            :src="asset('cspace_corner.png')"
            alt="Two abundance dimensions: field stars grey, cluster members magenta in one tight clump"
            style="width: 100%; height: auto"
          />
        </div>
      </div>
      <aside class="notes">
        (~2 min) Frame the input — this is the same dataset they will use in the hands-on. Make
        two points. First, standardisation is not cosmetic: every algorithm from here on is a
        statement about distance, and an unstandardised element with ten times the scatter would
        silently own that distance. Second, the plot is a trap — it shows 2 of 16 dimensions, and
        the red clump is only obvious because I coloured it. In the other 14 dimensions, and
        without the colours, nobody can see it. That gap is why the rest of the lecture exists.
      </aside>
    </section>

    <!-- 3 · The clustering problem -->
    <section>
      <div class="eyebrow">The problem</div>
      <h2>Finding groups in 16-D is not obvious</h2>
      <p class="small">
        In 2-D you can see an overdensity. In 16-D the volume explodes, "look for the peak" stops
        working, and clustering analysis earns its keep — under its own rules (Garcia-Dias et al. 2020).
      </p>
      <div class="cols" style="--n: 3; margin-top: 0.45em">
        <div class="panel">
          <h3>Multimodality first</h3>
          <p class="small">An algorithm <em>always</em> returns groups — even from a uniform cloud. Test for a multi-peaked distribution first (dip test).</p>
        </div>
        <div class="panel">
          <h3>Four tasks</h3>
          <p class="small">Feature selection → similarity metric → grouping criterion → validation. The grouping criterion is the core.</p>
        </div>
        <div class="panel flip">
          <h3>The trap</h3>
          <p class="small">"Distinct groups will be found even when there is no multipeak distribution" — the interpretation is on you.</p>
        </div>
      </div>
      <p class="small center muted" style="margin-top: 0.5em">
        There are no bad models — only models used outside their assumptions.
      </p>
      <aside class="notes">
        (~2 min) The chapter's core warning, and the BigDataLDN opening move. K-means (and most
        algorithms) will happily hand you clusters from pure noise. The honest question — for
        chemical tagging too — is not "did I get clusters?" but "are they real?" That's why
        validation gets its own section later. This slide also seeds the "four tasks": every
        algorithm that follows differs mainly in task 2 (metric) and task 3 (grouping criterion).
      </aside>
    </section>

    <!-- 4 · K-means — how it works -->
    <section>
      <div class="eyebrow">Tool one · K-means</div>
      <h2>K-means: the default tool</h2>
      <div class="fig-split" style="--cols: 1fr 1fr; margin-top: 0.3em; align-items: start">
        <div>
          <ol class="contribs small">
            <li><strong>Choose K</strong> and scatter K centres μ<sub>1</sub>…μ<sub>K</sub> at random.</li>
            <li><strong>Assign</strong> each star to its nearest centre (Euclidean).</li>
            <li><strong>Update</strong> each centre to the mean of its points.</li>
            <li><strong>Repeat</strong> 2–3 until no star changes cluster.</li>
          </ol>
          <p class="small" style="margin-top: 0.4em">
            Steps 2–3 descend a single objective, the within-cluster scatter
            <code>SSE = Σ<sub>j</sub> Σ<sub>i∈S<sub>j</sub></sub> ‖x<sub>i</sub> − μ<sub>j</sub>‖²</code>.
          </p>
          <p class="small muted" style="margin-top: 0.3em">
            Partitional · centroid-based · hard labels — every star ends up in exactly one
            cluster, and <strong>K is an input</strong>, never a result.
          </p>
        </div>
        <div class="figure">
          <img
            :src="asset('kmeans.gif')"
            alt="K-means iterating: points recoloured by nearest centroid, centroids moving each step"
            style="width: 100%; height: auto"
          />
        </div>
      </div>
      <aside class="notes">
        (~3 min) First algorithm of the day, so walk all four steps out loud while the GIF loops —
        point at the screen: recolour, move, recolour, move, stop. Say it is from the 1950s
        (Steinhaus 1956; Lloyd 1982; MacQueen 1967) and still the default first thing anyone
        reaches for. Then plant the two seeds that the next slides harvest: (a) K is an input, not
        an output — the algorithm cannot tell you how many clusters exist; (b) "nearest centre in
        Euclidean distance" means the boundaries between clusters are straight lines, so K-means
        can only carve the space into convex cells. Both come back as failure modes shortly.
      </aside>
    </section>

    <!-- 5 · Choose K and scatter K centres -->
    <section>
      <div class="eyebrow">Tool one · K-means · step 1 of 4</div>
      <h2>Choose K and scatter K centres</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('kmeans_step1_init.png')" alt="All points unassigned, K star-shaped centres dropped at random" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">K is an input: you decide how many clusters to look for, then drop K centres — here on K randomly chosen stars.</p>
      <aside class="notes">
        (~30 s) The only decision the algorithm cannot make for you. Land it: K is an input, never a
        result. Note the three centres already carry a colour and a shape — nothing is assigned
        yet, but each centre has an identity, so the colours on the next slide read as "belongs
        to that centre" rather than as three arbitrary groups.
      </aside>
    </section>

    <!-- 6 · Assign each star to its nearest centre -->
    <section>
      <div class="eyebrow">Tool one · K-means · step 2 of 4</div>
      <h2>Assign each star to its nearest centre</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('kmeans_step2_assign.png')" alt="Points coloured by their nearest centre; centres have not moved" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Each line ties one star to the centre that owns it — and where two bundles meet, the boundary is a straight line.</p>
      <aside class="notes">
        (~30 s) Trace one line with a finger: that star, that centre, nothing else. Then step back
        and let the bundles do the work. Two things to name. The seam where two fans meet is dead
        straight, because "nearest centre" is decided by a perpendicular bisector — that geometry
        is the failure mode we come back to. And these spokes are long, which is the point: add
        up their squared lengths and you have the SSE the algorithm is trying to shrink.
      </aside>
    </section>

    <!-- 7 · Move each centre to the mean of its stars -->
    <section>
      <div class="eyebrow">Tool one · K-means · step 3 of 4</div>
      <h2>Move each centre to the mean of its stars</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('kmeans_step3_update.png')" alt="Centres jumped to the mean of their assigned points; points have not moved" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Hollow marks show where each centre was; the arrow runs to the mean of the stars it owns. Not one star moved.</p>
      <aside class="notes">
        (~30 s) The update step, and the arrows are the whole point: from a bad start the centres
        travel a long way in one pass. Emphasise 'mean' — that is where the name comes from, and
        why a single outlier drags a centre. Say explicitly that the stars did not move and did
        not change colour; only the centres did.
      </aside>
    </section>

    <!-- 8 · Repeat until nothing changes -->
    <section>
      <div class="eyebrow">Tool one · K-means · step 4 of 4</div>
      <h2>Repeat until nothing changes</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('kmeans_step4_converged.png')" alt="Final stable partition: centres at cluster means, no point switching" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Eight passes later the spokes are as short as they can get — and their total squared length <em>is</em> the SSE.</p>
      <aside class="notes">
        (~30 s) Hold this next to step 2: there the spokes reach right across the field, here each
        centre sits at the heart of a short, tidy fan. That shortening is the algorithm working,
        because the total squared spoke length is exactly the SSE — the two moves are just the
        two ways of shortening it, and each can only make it smaller, so the loop must stop. It
        stops at a local minimum, which is the hook for the slide after next.
      </aside>
    </section>

    <!-- 9 · The loop, end to end -->
    <section>
      <div class="eyebrow">Tool one · K-means · the loop</div>
      <h2>Two moves, repeated until nothing changes</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 880px; margin: 0 auto">
        <img
          :src="asset('kmeans_storyboard.png')"
          alt="Six panels: the raw data, three centres dropped at random, the assignment drawn as spokes, the centres moving to their means, a second assignment, and the converged partition"
          style="width: 100%; height: auto"
        />
      </div>
      <p class="small muted center" style="margin-top: 0.3em; max-width: 880px; margin-inline: auto">(c)→(d) and (e)→(f) are the same two moves. That repetition <em>is</em> the algorithm.</p>
      <aside class="notes">
        (~30 s) The recap the four step slides cannot give, because you can never see two of
        them at once. Point at (c) and (e): same move. At (d) and (f): same move. Everything
        K-means does is those two, alternated, until the colours stop changing. Then set up the
        next slide by pointing back at (b): every one of these panels was decided by where those
        three centres happened to land.
      </aside>
    </section>

    <!-- 10 · K-means — initialization & SSE -->
    <section>
      <div class="eyebrow">Tool one · K-means</div>
      <h2>Run it again — you may get a worse answer</h2>
      <div class="fig-split" style="--cols: 1.05fr 0.95fr; margin-top: 0.3em; align-items: start">
        <div class="figure">
          <img
            :src="asset('kmeans_init.png')"
            alt="Two K-means runs from different initialisations on the same data: a good init finds the three-way split (SSE 331), a bad init merges two clusters (SSE 553)"
            style="width: 100%; height: auto"
          />
        </div>
        <div>
          <ul class="dotlist small">
            <li><strong>Initialisation matters.</strong> Same data, two seeds: the good run finds the three-way split (SSE = 331), the bad run merges two clusters (SSE = 553) — a worse <em>local</em> minimum.</li>
            <li><strong>SSE is not convex.</strong> K-means only reaches a local minimum: in the chapter, 4 of 5 runs score ≈89% against known labels, one collapses to ≈49%.</li>
            <li><strong>SSE</strong> (sum of squared error) discards those poor runs — but it barely separates the good ones, and it can never choose K.</li>
            <li><strong>Guardrails:</strong> K-means++ seeding (Arthur &amp; Vassilvitskii 2007) and many restarts — &gt;250 for one cortical segmentation (Nanetti et al. 2009).</li>
          </ul>
        </div>
      </div>
      <aside class="notes">
        (~3 min) Two ideas, and the audience must not conflate them. First, the figure: same
        data, two different seeds — the good run finds the three-way split, the bad run merges
        two clusters into one at a visibly higher SSE. That worse answer is a local minimum, not
        a bug (cluster labels still carry no meaning — never compare cluster "2" across runs).
        Second, the real hazard: SSE is a non-convex objective, so K-means descends into
        whichever local minimum its seed is nearest, and the chapter's example has one run in
        five landing at half the accuracy of the others. The cheap fix is restarts plus keeping
        the lowest SSE — K-means is fast enough that hundreds of restarts still beat one run of
        anything smarter. Close by flagging the gap this leaves: SSE cannot choose K, because it
        falls monotonically as K grows. Silhouette (Rousseeuw 1987) is the usual tool for that,
        and we come back to it in the validation section — with the same caveat as everything
        else here: only believe a clear peak.
      </aside>
    </section>

    <!-- 11 · K-means limits -->
    <section>
      <div class="eyebrow">Tool one · K-means limits</div>
      <h2>It assumes round, equal-sized clusters</h2>
      <div class="fig-split" style="--cols: 1.1fr 1fr; margin-top: 0.35em; align-items: start">
        <div class="figure" style="margin: 0">
          <img
            :src="asset('kmeans_fail.png')"
            alt="K-means splitting the two half-moons and cutting the big cluster in two"
            style="width: 100%; height: auto; display: block"
          />
        </div>
        <div>
          <ul class="crosslist small">
            <li><strong>Non-spherical shapes</strong> — a half-moon is split in two.</li>
            <li><strong>Different scales</strong> — one wide Gaussian eats the rest.</li>
            <li><strong>Unbalanced sizes</strong> — the small group is absorbed.</li>
            <li><strong>Straight boundaries</strong> — nearest centre cuts with lines.</li>
          </ul>
          <p class="small muted" style="margin-top: 0.4em">
            Give K-means the right K and it still fails here. Assigning each star to the nearest
            centre <em>in Euclidean distance</em> can only carve the space with straight edges. These
            are characteristics, not bugs — but they are exactly the geometry 16-D chemical space
            throws at us: uneven, anisotropic, overlapping groups.
          </p>
        </div>
      </div>
      <aside class="notes">
        (~2 min) Make it concrete with the chapter's own example (Fig. 13.4): give K-means the
        correct K and it still fails when the groups violate its assumptions — non-spherical
        shapes, different scales, unbalanced sizes, and the straight decision boundaries that come
        with Euclidean distance to a centre. Say the line out loud: these are characteristics, not
        bugs. That is the setup for the whole tour — every algorithm that follows exists to remove
        one of these assumptions.
      </aside>
    </section>

    <!-- 12 · KNN -->
    <section>
      <div class="eyebrow">Tool two · KNN</div>
      <h2>KNN: the nearest-neighbour primitive</h2>
      <p class="small">
        K-nearest neighbours asks one question — <em>which k points are closest to x?</em> It is not a
        clustering algorithm; it is the <strong>local primitive</strong> that DBSCAN, HDBSCAN*, UMAP and
        EVoC all call.
      </p>
      <div class="fig-split" style="--cols: 1.8fr 1fr; margin-top: 0.25em">
        <div class="panel">
          <h3>The algorithm</h3>
          <ol class="contribs tight small">
            <li><strong>Pick k</strong> — the neighbourhood size, and the method's only knob.</li>
            <li><strong>Measure</strong> the distance d(x, x<sub>i</sub>) to every other point.</li>
            <li><strong>Keep</strong> the k smallest → the k nearest neighbours of x.</li>
            <li><strong>Use them</strong> — a majority vote (classification), the k-th distance (density), or the edges themselves (graph).</li>
          </ol>
        </div>
        <div class="figure" style="margin: 0">
          <img
            :src="asset('knn_core_distance.png')"
            alt="Core distance: the k-th nearest-neighbour radius is small in dense regions and large in sparse ones"
            style="width: 100%; height: auto; display: block"
          />
        </div>
      </div>
      <div class="cols compact" style="--n: 2; margin-top: 0.25em">
        <div class="panel">
          <h3>Core distance κ(x)</h3>
          <p class="small">
            The distance to the k-th neighbour: small where stars are packed, large where they are
            sparse. Two of them give the <strong>mutual reachability</strong> distance
            max(κ(x), κ(y), d(x, y)) that HDBSCAN* clusters on.
          </p>
        </div>
        <div class="panel flip">
          <h3>The k-NN graph</h3>
          <p class="small">
            Keep the edges, not just the distance: every point joined to its k neighbours. That
            graph is the object <strong>UMAP</strong> and <strong>EVoC</strong> embed — clustering
            becomes a question about a graph.
          </p>
        </div>
      </div>
      <aside class="notes">
        (~4 min) The pivot slide — everything after this is a caller of it, so do not rush. KNN is
        normally taught as a supervised classifier, but here we steal its one useful idea: the
        neighbourhood. Two things fall out of it. First, the k-th-neighbour distance κ(x) is a free
        local density estimate — that is DBSCAN's minPts test and HDBSCAN*'s mutual reachability.
        Second, the neighbour lists are a graph, and that graph is what UMAP and EVoC embed. Point
        at the figure: same k, two very different radii, and that difference is the density signal.
      </aside>
    </section>

    <!-- 13 · DBSCAN -->
    <section>
      <div class="eyebrow">Tool three · DBSCAN</div>
      <h2>DBSCAN: follow the density</h2>
      <div class="slide-body">
      <div class="fig-split" style="--cols: 1.3fr 1fr; margin-top: 0.3em; align-items: start">
        <div>
          <p class="small">
            No centres at all: clusters are high-density regions, grown from neighbour counts
            instead of nearest-centre distance (Ester et al. 1996). Two parameters, three kinds of
            point.
          </p>
          <ol class="contribs tight small">
            <li><strong>Pick</strong> the neighbourhood radius <strong>ε</strong> and the minimum count <strong>minPts</strong>.</li>
            <li><strong>Core point:</strong> x is <em>core</em> if its ε-ball holds at least minPts points, counting x itself.</li>
            <li><strong>Grow:</strong> a cluster is a core point plus everything density-reachable from it — a chain of cores, each within ε of the next.</li>
            <li><strong>Border point:</strong> a non-core point inside a core's ε-ball joins that cluster, but never extends it.</li>
            <li><strong>Noise:</strong> reachable from no core → outlier. DBSCAN need not cluster everything.</li>
          </ol>
        </div>
        <div>
          <div class="figure" style="margin: 0">
            <img
              :src="asset('dbscan.gif')"
              alt="DBSCAN growing clusters: core points with radius epsilon, density-reachable points joining, outliers marked with a cross"
              style="width: 100%; height: auto; display: block"
            />
          </div>
          <p class="small muted" style="margin-top: 0.5em">
            The KNN thread: the core test is still a neighbour count — <strong>k fixed, radius
            free</strong> in KNN; <strong>radius fixed, k free</strong> in DBSCAN.
          </p>
        </div>
      </div>
      </div>
      <aside class="notes">
        (~3 min) Define both knobs out loud before anything else: ε is how far you look, minPts is
        how many you need to find. Then walk the animation — a core point recruits its ε-ball,
        the recruits that are themselves core keep the chain going, the ones that are not sit on
        the border and stop it. DBSCAN's two big wins over K-means: arbitrary shapes (no centroid,
        no Gaussian), and — uniquely so far — it labels outliers instead of forcing every point
        into a cluster. That "some points are noise" option is exactly right for chemical tagging,
        where most field stars belong to nothing. Close on the KNN thread: same neighbour count,
        only which of k and the radius is held fixed changes. But DBSCAN has its own assumption,
        coming next.
      </aside>
    </section>

    <!-- 14 · Pick ε and minPts -->
    <section>
      <div class="eyebrow">Tool three · DBSCAN · step 1 of 4</div>
      <h2>Pick ε and minPts</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('dbscan_step1_knobs.png')" alt="Two points, each with a circle of radius epsilon drawn around it: one holds five points and passes the minPts test, the other holds three and fails" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">ε is a radius, so draw it: put a ball of that size on a point and count what falls inside. Pass minPts and the point is <em>core</em>.</p>
      <aside class="notes">
        (~30 s) Two knobs, and both are visible here. ε is how far you look — literally the circle. minPts is how many you need to find inside it, counting the point itself. Show the pass and the near miss together: five clears the bar of four, three does not. Same neighbour-counting idea as KNN, only now the radius is fixed and the count varies, instead of the other way round.
      </aside>
    </section>

    <!-- 15 · Classify every point: core, border, noise -->
    <section>
      <div class="eyebrow">Tool three · DBSCAN · step 2 of 4</div>
      <h2>Classify every point: core, border, noise</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('dbscan_step2_classify.png')" alt="The same points, now labelled: twelve core points shown with their epsilon-balls, three border triangles inside a core's ball, four noise crosses inside nobody's" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Draw a ball on every core and the three roles read straight off: inside your own ball, inside someone else's, or inside nobody's.</p>
      <aside class="notes">
        (~30 s) One rule, applied to all nineteen points, and it sorts them into three piles. Core: the ball around it is crowded enough. Border: it fails the count itself, but it sits inside a core's ball, so it joins that cluster. Noise: no core's ball reaches it. Only the twelve cores can start or extend anything — that is the sentence to land before the next slide.
      </aside>
    </section>

    <!-- 16 · Grow clusters from core points -->
    <section>
      <div class="eyebrow">Tool three · DBSCAN · step 3 of 4</div>
      <h2>Link the cores, and the clusters appear</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('dbscan_step3_grow.png')" alt="Double-headed arrows between core points within epsilon of each other form two connected chains; single-headed arrows run out to the border triangles; noise points sit inside empty dashed balls" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Two cores within ε link <strong>both ways</strong>; a core links to a border <strong>one way only</strong>. A cluster is one connected chain — and nobody chose how many.</p>
      <aside class="notes">
        (~30 s) The arrows are the definition. Between two cores the link runs both ways, so the chain keeps going; out to a border it runs one way only, which is exactly why a border joins a cluster but can never grow it. Follow a chain and you have a cluster: two of them here, and nothing chose that number — it is however many connected chains the data happens to contain. That is the answer to K-means asking you for K.
      </aside>
    </section>

    <!-- 17 · Read off shapes and outliers -->
    <section>
      <div class="eyebrow">Tool three · DBSCAN · step 4 of 4</div>
      <h2>Read off shapes and outliers</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('dbscan_step4_result.png')" alt="DBSCAN on 300 points of interleaved half-moons plus uniform noise: both crescents recovered whole, with the scattered points marked as noise crosses" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Same two knobs, hundreds of stars — and these are the half-moons K-means cut down the middle on the limits slide.</p>
      <aside class="notes">
        (~30 s) Same ε, same minPts, no new machinery — just more points. Two wins over K-means, and both are on this one picture. Arbitrary shapes: these are the interleaved crescents that K-means sliced in half back on the limits slide, and a chain of ε-balls follows a curve as happily as a blob. And noise is an allowed answer — the scattered points are left out rather than forced into whichever cluster is nearest. For chemical tagging that second one is the point: most field stars belong to no cluster at all.
      </aside>
    </section>

    <!-- 18 · No free lunch -->
    <section>
      <div class="eyebrow">The toolbox, compared</div>
      <h2>No free lunch</h2>
      <p class="small">
        K-means and DBSCAN each win on their own toy data and lose on the others — a small change
        flips the ranking, and averaged over all problems no algorithm is universally best. KNN is
        not a clusterer: it is the <strong>primitive</strong> underneath DBSCAN and everything after.
      </p>
      <div class="cols" style="--n: 3; margin-top: 0.4em">
        <div class="panel">
          <h3>K-means</h3>
          <p class="small">Fast, simple. Assumes spherical, similar-size clusters; needs K; Euclidean → linear boundaries.</p>
        </div>
        <div class="panel flip">
          <h3>KNN</h3>
          <p class="small">The primitive. "Who's close?" — no clusters, just neighbours; feeds DBSCAN, HDBSCAN*, UMAP, EVoC.</p>
        </div>
        <div class="panel">
          <h3>DBSCAN</h3>
          <p class="small">Arbitrary shape + outliers. Assumes one density threshold works for every cluster.</p>
        </div>
      </div>
      <p class="stat center" style="margin-top: 0.4em">Know the assumptions</p>
      <p class="small center muted">The job is to match them to your data — and to check with validation.</p>
      <aside class="notes">
        (~2 min) The no-free-lunch beat (the BigDataLDN through-line). No model is "bad" — each
        encodes assumptions, and the failure is applying a model outside them. KNN is the odd one
        out by design: it makes no clusters at all, it just measures locality — which is why it
        shows up inside every method that follows. This is the philosophical core the whole school
        wants you to internalise, and it leads straight into the density section.
      </aside>
    </section>

    <!-- 19 · Validation — internal -->
    <section>
      <div class="eyebrow">Validation · internal scores</div>
      <h2>How do you know the clusters are real?</h2>
      <p class="small">
        <strong>Internal</strong> validation uses only the data and the labels the algorithm just
        produced — no answer key. Three scores, three different questions, and a sharp caveat on
        each.
      </p>
      <div class="cols" style="--n: 3; margin-top: 0.4em">
        <div class="panel">
          <h3>SSE</h3>
          <p class="small">
            Summed squared distance from every object to its own cluster mean. It falls with every
            extra cluster, so it cannot choose K — but across random restarts the run with the
            lowest SSE is the one to keep.
          </p>
        </div>
        <div class="panel">
          <h3>Silhouette</h3>
          <p class="small">
            Mean within-cluster distance against mean between-cluster distance, averaged over every
            object: −1 to 1, where 1 is cleanly separated and 0 is completely overlapping. Its peak
            over K suggests K — trust it only when the peak is clear.
          </p>
        </div>
        <div class="panel flip">
          <h3>Dip test</h3>
          <p class="small">
            Run this <em>first</em>: is the distribution multipeaked at all? The null hypothesis is a
            single peak, so a small p-value means real structure. Silhouette only means something
            once this passes (Hartigan &amp; Hartigan 1985).
          </p>
        </div>
      </div>
      <p class="small center muted" style="margin-top: 0.5em">
        Homogeneity and ARI need true labels — that is <em>external</em> validation, next slide.
      </p>
      <aside class="notes">
        (~3 min) The chapter's §13.2.2, in the order you should actually use them: (1) is the data
        multipeaked at all — the dip test, and if it fails, stop; (2) how many groups — the
        silhouette, and only if it has a clear peak, because on unstructured data it returns
        essentially random values; (3) which run to keep — the lowest SSE across restarts, since
        SSE falls monotonically with K and can never pick K for you. Say the honest bottom line the
        chapter states: these scores are aids, not oracles, and the interpretation has to be
        informed by domain knowledge. Then flag the trap on the last line — homogeneity (Rosenberg
        &amp; Hirschberg 2007, and the score the chapter itself uses to compare K-means, GMM and
        DBSCAN) and ARI get quoted as if they were internal scores, but both compare against truth,
        which is a completely different kind of check. The luxury of this problem is that we have
        one, which is the next slide.
      </aside>
    </section>

    <!-- 20 · Validation — external truth -->
    <section>
      <div class="eyebrow">Validation · external truth</div>
      <h2>The luxury of this problem: ground truth</h2>
      <p class="small">
        Internal scores can only say a grouping is tight and separated. <strong>External</strong>
        validation asks the harder question — is it the <em>right</em> grouping? — and that needs
        labels the clustering never saw.
      </p>
      <div class="cols" style="--n: 3; margin-top: 0.4em">
        <div class="panel">
          <h3>Where it comes from</h3>
          <p class="small">
            Chemistry is not the only fingerprint: stars born together also share position,
            parallax, proper motion and radial velocity. Gaia kinematics label membership without
            ever touching the 16 abundances we clustered on.
          </p>
        </div>
        <div class="panel">
          <h3>What it buys</h3>
          <p class="small">
            <strong>Recall</strong> — did we recover the real members? <strong>Precision</strong> — are
            the ones we claim real? And homogeneity and ARI, the label-invariant scores, finally
            mean something.
          </p>
        </div>
        <div class="panel flip">
          <h3>Why it matters</h3>
          <p class="small">
            Most unsupervised problems — brain imaging, say — have no answer key, so the argument
            never ends. Chemical tagging has one, and that is what turns it into a benchmark: the
            results later in this talk are measured, not hoped for.
          </p>
        </div>
      </div>
      <p class="small center muted" style="margin-top: 0.5em">
        Internal asks <em>is this grouping self-consistent?</em> — external asks <em>is this grouping true?</em>
      </p>
      <aside class="notes">
        (~2 min) The pivot from "trust the score" to "check against reality", and the contrast with
        the previous slide is the point: internal scores never leave the data you clustered,
        external validation brings in an answer from outside it. Stress the independence — the
        labels come from six-dimensional phase space, not from the 16 abundances, so scoring
        chemistry against kinematics is not circular. That is why the benchmark at the end means
        something: we do not have to argue about whose silhouette is better, we can quote recall
        and precision. Keep this in your back pocket — it returns at the benchmark slides.
      </aside>
    </section>

    <!-- 21 · The scale problem -->
    <section>
      <div class="eyebrow">The scale problem</div>
      <h2>One density threshold can't see everything</h2>
      <p class="small">
        DBSCAN commits to a single density threshold ε. But stellar structure lives across a
        <strong>huge range of densities</strong>, and one ε can only ever pick one point on that axis.
      </p>
      <p class="small muted" style="margin: 0.45em 0 0.1em">Sparse → dense, all in the same survey:</p>
      <div class="scale-bar"></div>
      <div class="scale-ticks">
        <span>tidal tails</span>
        <span>moving groups</span>
        <span>the field disc</span>
        <span>open clusters</span>
        <span>globular cores</span>
      </div>
      <div class="cols" style="--n: 3; margin-top: 0.55em">
        <div class="panel">
          <h3>The symptom</h3>
          <p class="small">
            Loosen ε until the tidal tail hangs together, and the disc fuses into one blob. Tighten it
            until the globular splits cleanly, and the tail is all noise.
          </p>
        </div>
        <div class="panel">
          <h3>The half-fix</h3>
          <p class="small">
            <strong>OPTICS</strong> (Ankerst et al. 1999) refuses to choose: it plots the density profile at every
            scale and lets the analyst read the valleys. Honest — but still a picture a human has to
            interpret, one dataset at a time.
          </p>
        </div>
        <div class="panel flip">
          <h3>What we want</h3>
          <p class="small">
            A <em>hierarchy over all density scales</em>, with every cluster extracted at the scale where it
            actually exists — and the reading-off done for us.
          </p>
        </div>
      </div>
      <aside class="notes">
        (~2 min) The transition slide. Say it plainly: DBSCAN's one assumption — a single global
        density threshold — is the first thing real data breaks. Walk the bar left to right and name
        a real object at each end so they feel the dynamic range: a dissolving tidal tail is orders
        of magnitude sparser than a globular core, and no single ε serves both. Then OPTICS as the
        honest halfway house — it shows you the whole profile instead of picking for you, but a
        human still has to read the picture, and that does not scale to a survey. The next two
        slides automate exactly that reading: HDBSCAN* builds the hierarchy, PLSCAN selects from it.
      </aside>
    </section>

    <!-- 22 · HDBSCAN* -->
    <section>
      <div class="eyebrow">Density at all scales · HDBSCAN*</div>
      <h2>HDBSCAN*: a hierarchy over every density</h2>
      <p class="small">
        HDBSCAN* (Campello et al. 2013, 2015; McInnes &amp; Healy 2017) runs DBSCAN at
        <strong>every threshold at once</strong> and keeps the whole tree of answers.
      </p>
      <div class="cols compact" style="--n: 2; margin-top: 0.25em">
        <div class="panel">
          <h3>The machinery</h3>
          <p class="small">
            <strong>Mutual reachability.</strong> d<sub>mut</sub> = max(κ<sub>i</sub>, κ<sub>j</sub>, d<sub>ij</sub>), with κ the distance to the
            k-th neighbour. Gaps shorter than κ are levelled up, pushing sparse points away.
          </p>
          <p class="small">
            <strong>Build, then read.</strong> Single-linkage on d<sub>mut</sub> gives one branch per structure over
            the density range it survives; slicing at one density recovers DBSCAN.
          </p>
        </div>
        <div class="panel flip">
          <h3>The knob you still pick</h3>
          <p class="small">
            Two survive: <strong>k</strong> (behind κ) and <strong>min_cluster_size m<sub>c</sub></strong> — how many stars a branch
            needs to count as a cluster rather than a wiggle.
          </p>
          <p class="small">
            m<sub>c</sub> is a <strong>smoothing knob</strong>: raise it and shallow peaks disappear, as if the density
            were blurred (Bot et al. 2025, Fig. 1). m<sub>c</sub> = 5 makes a five-star clump a cluster; m<sub>c</sub> =
            50 makes it noise.
          </p>
        </div>
      </div>
      <div class="fig-split" style="--cols: 1.54fr 1.29fr; max-width: 640px; margin: 0.25em auto 0">
        <div class="figure" style="width: 100%; box-sizing: border-box">
          <img :src="asset('mutual_reachability.png')" alt="Two points, each ringed by its core-distance circle, with the straight-line distance between them and the max-of-three definition of mutual reachability" style="width: 100%; height: auto" />
        </div>
        <div class="figure" style="width: 100%; box-sizing: border-box">
          <img :src="asset('hdbscan_density.gif')" alt="Animation sweeping the density level: clusters appear, grow and merge as the threshold drops" style="width: 100%; height: auto" />
        </div>
      </div>
      <aside class="notes">
        (~3 min) Two ideas, in this order. First mutual reachability, pointing at the left figure:
        taking the max of the two core distances and the true distance means any gap smaller than a
        point's own core distance gets levelled up to it. Sparse points end up far from everything,
        which is exactly what stops single-linkage chaining through noise. Second, the sweep on the
        right — narrate it as "this is DBSCAN, at every ε, simultaneously". The tree records the
        level at which each clump is born and the level at which it merges away, and the branch that
        survives the longest span is the one you want. Then set up the next slide: the density
        threshold is gone, but m_c is not, and m_c is a genuine judgement call with no data-driven
        answer. (The right-hand panel is an animation — if you are presenting from the PDF export,
        talk through the sweep instead of pointing at it.)
      </aside>
    </section>

    <!-- 23 · Core distance κ(x) -->
    <section>
      <div class="eyebrow">Tool four · HDBSCAN* · step 1 of 4</div>
      <h2>Core distance κ(x)</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 740px; margin: 0 auto">
        <img :src="asset('knn_core_distance.png')" alt="Two query points, each ringed by its core-distance circle" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">The distance to the k-th neighbour — the local density probe you already met in KNN.</p>
      <aside class="notes">
        (~30 s) Reuse the KNN slide's picture deliberately: κ is the same quantity. Dense regions give small κ, sparse regions large κ.
      </aside>
    </section>

    <!-- 24 · Mutual reachability d_mut -->
    <section>
      <div class="eyebrow">Tool four · HDBSCAN* · step 2 of 4</div>
      <h2>Mutual reachability d_mut</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('mutual_reachability.png')" alt="Two points with core-distance circles and the max-of-three definition" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">d_mut = max(κ_i, κ_j, d_ij): gaps shorter than κ are levelled up, sparse points are pushed apart.</p>
      <aside class="notes">
        (~30 s) The key move. Taking the max flattens dense regions and separates sparse points — this is what stops single-linkage from chaining through noise.
      </aside>
    </section>

    <!-- 25 · Build one tree over all densities -->
    <section>
      <div class="eyebrow">Tool four · HDBSCAN* · step 3 of 4</div>
      <h2>Build one tree over all densities</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('hdbscan_step3_tree.png')" alt="A single-linkage dendrogram with a horizontal density cut line" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Single-linkage on d_mut gives one branch per structure — slicing at any density recovers DBSCAN.</p>
      <aside class="notes">
        (~30 s) One tree contains every DBSCAN run. The horizontal cut is one density threshold; the tree keeps all of them.
      </aside>
    </section>

    <!-- 26 · Keep the longest-lived branches -->
    <section>
      <div class="eyebrow">Tool four · HDBSCAN* · step 4 of 4</div>
      <h2>Keep the longest-lived branches</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('hdbscan_step4_sweep.png')" alt="Clusters visible at one density level of the sweep" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Run every density threshold at once; the branches that survive the longest span are your clusters.</p>
      <aside class="notes">
        (~30 s) The sweep is DBSCAN at every ε simultaneously. Longest-lived branch wins — that is the only selection rule, and m_c is the knob left over for next slide.
      </aside>
    </section>

    <!-- 27 · PLSCAN -->
    <section>
      <div class="eyebrow">Persistence · PLSCAN</div>
      <h2>PLSCAN: drop the min-cluster-size knob</h2>
      <p class="small">
        <strong>PLSCAN</strong> — Persistent Leaves Spatial Clustering for Applications with Noise (Bot,
        McInnes &amp; Aerts 2025) — keeps HDBSCAN*'s hierarchy but stops picking a scale. Instead it
        measures <strong>how long each cluster survives across all of them</strong>.
      </p>
      <div class="cols compact" style="--n: 3; margin-top: 0.45em">
        <div class="panel">
          <h3>1 · Scale-space</h3>
          <p class="small">
            Raising m<sub>c</sub> never moves a merge — it only prunes branches too small to count. So one
            condensed tree already contains every HDBSCAN* run at every m<sub>c</sub>.
          </p>
        </div>
        <div class="panel">
          <h3>2 · Leaf tree</h3>
          <p class="small">
            For each leaf cluster, record the interval of m<sub>c</sub> over which it is a leaf. Its length,
            s<sub>max</sub> − s<sub>min</sub>, is that cluster's <strong>persistence</strong> — one bar of the barcode.
          </p>
        </div>
        <div class="panel flip">
          <h3>3 · Persistence trace</h3>
          <p class="small">
            Sum the persistence of every leaf alive at each m<sub>c</sub>. Peaks in that trace are the most
            stable scales, and PLSCAN returns one <strong>layer of detail</strong> per peak.
          </p>
        </div>
      </div>
      <p class="small muted" style="margin: 0.45em 0 0">
        Formally, 0-dimensional persistent homology on a purpose-built metric space (their App. A).
        Practically: no m<sub>c</sub> to guess. <strong>k survives</strong> — which is the next slide's business.
      </p>
      <div class="fig-split" style="--cols: 1.42fr 2.39fr; max-width: 900px; margin: 0.4em auto 0">
        <div class="figure" style="width: 100%; box-sizing: border-box">
          <img :src="asset('plscan_persistence.gif')" alt="Sweeping min-cluster-size across a three-peak density profile; each peak's bar grows for as long as that cluster stays alive" style="width: 100%; height: auto" />
        </div>
        <div class="figure" style="width: 100%; box-sizing: border-box">
          <img :src="asset('plscan_barcode.png')" alt="A three-peak density profile beside its persistence barcode: one bar per cluster, its length the range of min-cluster-size over which that cluster exists" style="width: 100%; height: auto" />
        </div>
      </div>
      <aside class="notes">
        (~3 min) The newest material in the talk, so go slowly and lean on the pictures. Start from
        last slide's complaint: m_c is a smoothing knob with no right answer. PLSCAN's move is to
        refuse the question — sweep m_c over its whole range and ask which clumps are still there at
        the end. Run the left animation and watch peaks drop out as the threshold rises. The right
        figure is the same information redrawn: one bar per cluster, bar length = the span of m_c
        over which it survived = its persistence. Here the densest peak persists to 202 and the
        shallowest only to 70 — say that a bar this short would be discarded as a fluctuation on
        real data, but note the toy figure only draws the three survivors, not the stubs. Be honest
        about the caveat in the muted line: PLSCAN is not parameter-free, k is still there. What it
        removes is the parameter you had no principled way to set.
      </aside>
    </section>

    <!-- 28 · Find the density peaks -->
    <section>
      <div class="eyebrow">Tool five · PLSCAN · step 1 of 3</div>
      <h2>Find the density peaks</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('plscan_step1_density.png')" alt="A three-peak density profile with each peak marked" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Every density peak is a candidate cluster — no min_cluster_size to guess.</p>
      <aside class="notes">
        (~30 s) The starting point: peaks in the density profile. PLSCAN refuses to pick a scale and instead considers all of them.
      </aside>
    </section>

    <!-- 29 · Measure each cluster's persistence -->
    <section>
      <div class="eyebrow">Tool five · PLSCAN · step 2 of 3</div>
      <h2>Measure each cluster's persistence</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('plscan_step2_barcode.png')" alt="One horizontal bar per cluster, length equal to its persistence" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">One bar per cluster: its length is how long that cluster survives as min_cluster_size sweeps.</p>
      <aside class="notes">
        (~30 s) Persistence = lifetime in min_cluster_size. A bar this long means the cluster is still there across a wide range of scales.
      </aside>
    </section>

    <!-- 30 · Keep the long bars, discard the stubs -->
    <section>
      <div class="eyebrow">Tool five · PLSCAN · step 3 of 3</div>
      <h2>Keep the long bars, discard the stubs</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('plscan_step3_read.png')" alt="Barcode with long bars coloured and short bars greyed, a cut line between" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Long bars are stable clusters; short bars are fluctuations — the data, not you, draws the line.</p>
      <aside class="notes">
        (~30 s) The long/short split replaces the m_c knob you used to set by hand. Caveat: k still survives — PLSCAN is not parameter-free.
      </aside>
    </section>

    <!-- 31 · PLSCAN — results & segue -->
    <section>
      <div class="eyebrow">Persistence · PLSCAN</div>
      <h2>More stable, less sensitive — and a familiar name</h2>
      <div class="slide-body">
      <p class="small">
        Bot et al. benchmark PLSCAN against HDBSCAN* (excess-of-mass selection) across a suite of
        real-world datasets, at the conventional default k = 4. Medians over those datasets:
      </p>
      <div class="kpis" style="--n: 4; margin-top: 0.45em">
        <div class="panel center">
          <div class="stat">0.66</div>
          <div class="stat-label">median ARI<br />HDBSCAN*: 0.57</div>
        </div>
        <div class="panel center">
          <div class="stat">0.85</div>
          <div class="stat-label">V-measure<br />HDBSCAN*: 0.79</div>
        </div>
        <div class="panel center">
          <div class="stat pink">0.74</div>
          <div class="stat-label">homogeneity<br />HDBSCAN*: 0.82</div>
        </div>
        <div class="panel center">
          <div class="stat pink">0.76</div>
          <div class="stat-label">non-noise fraction<br />HDBSCAN*: 0.88</div>
        </div>
      </div>
      <div class="cols compact" style="--n: 3; margin-top: 0.5em">
        <div class="panel">
          <h3>Less sensitive to k</h3>
          <p class="small">
            HDBSCAN*'s labels swing about at low k. PLSCAN returns much the same clustering across
            the tested range.
          </p>
        </div>
        <div class="panel">
          <h3>Affordable</h3>
          <p class="small">
            Run-times competitive with K-means++ on low-dimensional data; at high dimension it scales
            like HDBSCAN*, which is the price of using a space tree at all.
          </p>
        </div>
        <div class="panel flip">
          <h3>The trade-off</h3>
          <p class="small">
            The pink numbers are the bill: PLSCAN's clusters are more complete but less pure, and it
            sends more stars to the noise bin. A better ARI is not "better everywhere".
          </p>
        </div>
      </div>
      <p class="small center muted" style="margin-top: 0.5em">
        Middle author: <strong>Leland McInnes</strong>, also behind UMAP and the HDBSCAN* implementation
        everyone actually runs. Which is where we go next.
      </p>
      </div>
      <aside class="notes">
        (~2 min) Read the four numbers as one sentence, not four: PLSCAN wins on the summary scores
        because it recovers more of each true cluster, and it pays for that with lower purity and
        more stars thrown to noise. That is the deck's through-line again — no free lunch, only
        different assumptions. Do not oversell it: this is one benchmark suite at one default k, on
        general-purpose datasets rather than on abundances, and the authors say themselves that
        whether density maxima are the "true" clusters depends on the use case. Then the segue: the
        middle author is the person behind UMAP and behind the HDBSCAN* library, which is why EVoC
        later fuses all three — one research programme, not three coincidences. And now we change
        strategy entirely: stop clustering in the native space, and reshape the space first.
      </aside>
    </section>

    <!-- 32 · t-SNE — project to 2-D -->
    <section>
      <div class="eyebrow">Embeddings · t-SNE</div>
      <h2>t-SNE: project to 2-D, then look</h2>
      <p class="small">
        Everything so far clustered in the native abundance space. Kos et al. (2017) change the
        question: <strong>reshape the space first</strong>. We are excellent at spotting a clump in ≤ 3-D and
        hopeless in 13-D, so t-distributed stochastic neighbour embedding compresses the chemistry
        into a picture — and clustering becomes drawing a line round what you can already see.
      </p>
      <div class="fig-split" style="--cols: 1fr 1fr; align-items: start; margin-top: 0.5em">
        <div class="panel">
          <h3>Their recipe, on GALAH</h3>
          <p class="small">
            <strong>13 abundances, Manhattan distance.</strong> Summing |Δ| instead of squaring it, so one bad
            line does less damage than under Euclidean.
          </p>
          <p class="small">
            <strong>Weights ∝ 1 / cluster scatter.</strong> Ba and K, whose scatter is mostly noise, drop to
            0.25; the clean Fe, Ti, Cr and Cu rise to 2.0 — a stand-in for per-star errors they never
            measured.
          </p>
          <p class="small">
            <strong>Perplexity.</strong> How many neighbours count as "local" — next slide.
          </p>
          <p class="small">
            <strong>One sky region per cluster, 30–45° in radius.</strong> 9,408 stars within 40° of the
            Pleiades, say — never the whole survey at once.
          </p>
        </div>
        <div class="figure" style="width: 100%; box-sizing: border-box">
          <img :src="asset('tsne.gif')" alt="A t-SNE embedding converging from a scattered layout into four separated clusters" style="width: 100%; height: auto" />
        </div>
      </div>
      <aside class="notes">
        (~3 min) Flag the change of strategy first — it is the hinge of the lecture. Up to now we
        asked an algorithm to find clusters in high dimension; from here we compress to two and hand
        the problem to the audience's visual cortex. Then the three practical choices that are easy
        to skim past and matter enormously. Manhattan distance for outlier robustness. The weights,
        which are the inverse of the observed cluster scatter — Ba and K scatter mostly because they
        are hard to measure, so they get down-weighted, while Fe, Ti, Cr and Cu are trusted; it is a
        poor man's error bar, and Kos et al. say the map falls apart without it. And above all the
        region cut: they never ran one all-sky t-SNE, every map is a 30–45° cone around one known
        cluster. Promise them the region cut comes back to bite us in the benchmark.
      </aside>
    </section>

    <!-- 33 · t-SNE — how it works -->
    <section>
      <div class="eyebrow">Embeddings · t-SNE</div>
      <h2>t-SNE, under the hood</h2>
      <p class="small">
        t-SNE never tries to preserve distances. It turns both spaces into probability distributions
        over "who is whose neighbour", then matches them.
      </p>
      <div class="cols compact" style="--n: 3; grid-template-columns: 1.3fr 1fr 1fr; margin-top: 0.5em">
        <div class="panel">
          <h3>The mechanics</h3>
          <ul class="dotlist small">
            <li>High-D: Gaussian p<sub>ij</sub>, one bandwidth σ<sub>i</sub> per star</li>
            <li>Low-D: heavy-tailed Student-t q<sub>ij</sub> on the map</li>
            <li>Minimise <strong>KL(P‖Q)</strong> by gradient descent</li>
            <li>Barnes-Hut brings the cost to O(N log N)</li>
            <li>Random start, so every run differs</li>
            <li>Kos et al. keep the lowest-KL map of many</li>
          </ul>
        </div>
        <div class="panel">
          <h3>Perplexity: the knob</h3>
          <p class="small">
            Each σ<sub>i</sub> is tuned until star i's neighbour distribution has the <strong>perplexity</strong> you
            asked for — an effective neighbour count, typically 5–50. It plays the role k plays in
            KNN.
          </p>
          <p class="small">
            Low → fine local clumps. High → coarse global shape. Same data, different picture, so
            always quote it with the map.
          </p>
        </div>
        <div class="panel flip">
          <h3>The crowding problem</h3>
          <p class="small">
            A 13-D ball has far more room at middling distances than any 2-D disc. Match Gaussians in
            both spaces and all those middling neighbours pile inwards into one blob.
          </p>
          <p class="small">
            The Student-t's <strong>heavy tail</strong> is the fix: a decent q<sub>ij</sub> is reachable at a long map
            distance, so points can spread and gaps open between clumps.
          </p>
        </div>
      </div>
      <aside class="notes">
        (~2 min) The one mechanics slide; keep it to two takeaways. Takeaway one is perplexity: it is
        the local/global dial, roughly "how many neighbours count as near", and a map is meaningless
        without it quoted. Give the 5–50 range and point people at the Wattenberg demo that Kos et
        al. cite. Takeaway two is the crowding problem, because it explains the "t" in the name: if
        the map used a Gaussian too, there simply is not enough room in two dimensions for
        everything that sits at moderate distance in thirteen, and the whole plot collapses inward.
        The heavy tail buys that room back, and the visible gaps between clumps are its doing.
        Mention the random initialisation in passing — it is why two runs look different, and why
        anyone showing you a t-SNE map owes you the perplexity. Save the harder criticisms (cluster
        sizes and inter-cluster distances mean nothing) for the limits slide.
      </aside>
    </section>

    <!-- 34 · Model who is whose neighbour, in high-D -->
    <section>
      <div class="eyebrow">Embeddings · t-SNE · step 1 of 3</div>
      <h2>Model who is whose neighbour, in high-D</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('tsne_step1_highd.png')" alt="Two Gaussian similarity curves, narrow and wide" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Each star gets a Gaussian over its neighbours, tuned to a fixed perplexity — an effective neighbour count.</p>
      <aside class="notes">
        (~30 s) High-D similarities p_ij. Perplexity plays the role k played in KNN: how many neighbours count as near.
      </aside>
    </section>

    <!-- 35 · Model the same in 2-D, with a heavy tail -->
    <section>
      <div class="eyebrow">Embeddings · t-SNE · step 2 of 3</div>
      <h2>Model the same in 2-D, with a heavy tail</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('tsne_step2_lowd.png')" alt="A Gaussian curve beside a heavy-tailed Student-t curve" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">A Student-t in the map — the heavy tail lets points spread instead of crowding into one blob.</p>
      <aside class="notes">
        (~30 s) The crowding problem and its fix. The 't' in t-SNE is this Student-t; the visible gaps between clumps are its doing.
      </aside>
    </section>

    <!-- 36 · Match the two, by gradient descent -->
    <section>
      <div class="eyebrow">Embeddings · t-SNE · step 3 of 3</div>
      <h2>Match the two, by gradient descent</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('tsne_step3_descent.png')" alt="Final t-SNE embedding with four separated clusters" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Minimise KL(P‖Q): similar stars are pulled together, dissimilar ones pushed apart.</p>
      <aside class="notes">
        (~30 s) Gradient descent matches the two distributions. Random start, so every run differs — quote the perplexity with the map.
      </aside>
    </section>

    <!-- 37 · t-SNE results -->
    <section>
      <div class="eyebrow">t-SNE results · Kos et al. 2017</div>
      <h2>It recovered clusters — and found new Pleiades</h2>
      <p class="small muted" style="margin: 0.1em 0 0.3em">
        GALAH + K2-HERMES, <strong>13 abundances</strong>, 9408 stars in a 40° radius around the
        Pleiades. t-SNE draws the map; the groups are then drawn on it <em>by eye</em>.
      </p>
      <div class="fig-split" style="--cols: 1fr 1.4fr; align-items: start">
        <div class="figure-placeholder" style="aspect-ratio: 4 / 3; min-height: 0">
          <p class="figure-placeholder__tag">figure · slot in</p>
          <p class="figure-placeholder__desc">t-SNE map of the Pleiades region: members red, field grey, the two new members marked (paper's Fig. 2).</p>
          <p class="figure-placeholder__target">save as <code>kos_pleiades_tsne.png</code></p>
        </div>
        <div>
          <ul class="checklist small">
            <li><strong>7 of 9</strong> clusters recovered from chemistry alone</li>
            <li>Nine observed: six globular, three open</li>
            <li><strong>Two new Pleiades members</strong>, kinematically confirmed</li>
            <li>One of them 6° out — a full tidal radius</li>
            <li>Sub-groups of a cluster land together on the map</li>
          </ul>
          <ul class="crosslist small">
            <li><strong>47 Tuc</strong> — no chemical group at all; too like the field</li>
            <li><strong>NGC 2516</strong> — only an edge observed, so no verdict</li>
            <li>Field-star contamination inside every group</li>
            <li><em>“a large fraction of the stars are untaggable”</em></li>
          </ul>
        </div>
      </div>
      <aside class="notes">
        (~3 min) This is the proof the approach works on real stars — say the numbers out loud, then
        spend equal time on the right-hand column. Two of the nine failed for completely different
        reasons: 47 Tuc simply has no chemical group in the map, while NGC 2516 was barely observed,
        so it is a data failure, not a method failure. Worth flagging that the Pleiades sub-groups A
        and B the paper sees in The Cannon abundances could not be confirmed with SME — the
        hierarchy claim is about the projection preserving structure, not about that particular
        split being real. Close on the paper's own verdict: with 13 elements a large fraction of
        stars are untaggable. Chemical tagging is hard — which is exactly why we benchmark it.
      </aside>
    </section>

    <!-- 38 · t-SNE limits -->
    <section>
      <div class="eyebrow">t-SNE limits</div>
      <h2>t-SNE embeds — it does not cluster</h2>
      <div class="cols" style="--n: 2; margin-top: 0.3em">
        <div class="panel">
          <h3>What it won't tell you</h3>
          <ul class="crosslist small">
            <li><strong>Distance</strong> — gaps between clumps are arbitrary</li>
            <li><strong>Density</strong> — area on the map isn't volume</li>
            <li><strong>Size</strong> — sparse groups spread, dense ones shrink</li>
            <li><strong>Stability</strong> — rerun it and the clumps move</li>
            <li><strong>Speed</strong> — O(N²), even with Barnes-Hut</li>
          </ul>
        </div>
        <div class="panel flip">
          <h3>Where do labels come from?</h3>
          <p class="small">
            Not from t-SNE. It hands you a <em>picture</em> — Kos et al. drew the polygons by hand.
            To automate that you bolt a clusterer onto the map: <strong>HDBSCAN* on the 2-D
            embedding</strong>.
          </p>
          <ul class="crosslist small">
            <li>Two algorithms, two sets of knobs</li>
            <li>It clusters the picture, not the data</li>
            <li>So it inherits every distortion on the left</li>
          </ul>
        </div>
      </div>
      <p class="small center" style="margin-top: 0.4em">
        Read a t-SNE map for <strong>neighbourhoods</strong> — never distance, size or density.
      </p>
      <p class="small center muted" style="margin-top: 0.15em">
        Can one method embed <em>and</em> cluster — and pick the scale itself?
      </p>
      <aside class="notes">
        (~2 min) This is the slide to slow down on: the single most common mistake with t-SNE is
        reading it as a map with a scale. Say it plainly — t-SNE embeds, it does not cluster, and
        neither the distances nor the densities in the picture are trustworthy. The mitigation is
        the two-step pipeline, which works but inherits every weakness of both halves, because
        HDBSCAN* is clustering the distortion, not the abundances. The last line sets up the two
        requirements EVoC will eventually meet: one method, and a scale it chooses itself.
      </aside>
    </section>

    <!-- 39 · UMAP -->
    <section>
      <div class="eyebrow">Embeddings · UMAP</div>
      <h2>UMAP: a graph embedding with global structure</h2>
      <p class="small muted" style="margin: 0.1em 0 0.3em">
        Uniform Manifold Approximation and Projection (McInnes, Healy &amp; Melville 2018) skips the
        N² pairs entirely: it builds a <strong>weighted neighbour graph</strong>, then draws that
        graph. Hold on to the graph — EVoC reuses it, and so does its author's other work.
      </p>
      <div class="fig-split" style="--cols: 1.6fr 1fr; align-items: start">
        <div>
          <ol class="contribs small tight">
            <li><strong>kNN graph</strong> — join each star to its <strong>n_neighbors</strong> neighbours.</li>
            <li><strong>Fuzzy edges</strong> — rescale each point's distances by its own nearest-neighbour distance, so a weight is the <em>probability an edge exists</em>, not a distance.</li>
            <li><strong>Layout</strong> — attract along edges, repel sampled non-neighbours, minimising cross-entropy against that graph (t-SNE minimises KL).</li>
          </ol>
          <ul class="dotlist small" style="margin-top: 0.35em">
            <li><strong>n_neighbors</strong> — how local the graph is (cf. perplexity)</li>
            <li><strong>min_dist</strong> — how tightly the layout may pack points</li>
          </ul>
        </div>
        <div>
          <div class="figure">
            <img :src="asset('umap.gif')" alt="UMAP embedding converging from a spectral layout to four separated clusters" style="width: 100%; height: auto" />
          </div>
          <ul class="checklist small" style="margin-top: 0.4em">
            <li>Near-linear in N</li>
            <li>Keeps global layout</li>
            <li>Still a projector</li>
            <li>HDBSCAN* still after</li>
          </ul>
        </div>
      </div>
      <aside class="notes">
        (~3 min) Walk the three numbered steps slowly — this is the one place in the talk where a
        graph, not a distance matrix, is the object being optimised, and everything after it depends
        on that. The fuzzy step is the one to dwell on: because each point's distances are rescaled
        by its own nearest neighbour, a sparse region and a dense region get comparable edge
        weights, which is how UMAP avoids t-SNE's density distortion. Then the payoff: it is
        near-linear rather than quadratic, and unlike t-SNE the layout of clusters relative to each
        other means something. But it is still only a projector — HDBSCAN* still has to run
        afterwards. Land the last line hard: EVoC starts from this same graph.
      </aside>
    </section>

    <!-- 40 · Build the k-nearest-neighbour graph -->
    <section>
      <div class="eyebrow">Embeddings · UMAP · step 1 of 3</div>
      <h2>Build the k-nearest-neighbour graph</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('umap_step1_graph.png')" alt="Points joined by edges to their nearest neighbours" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Join each star to its k neighbours — the graph is the object everything after this is drawn on.</p>
      <aside class="notes">
        (~30 s) The same graph EVoC will reuse. Unlike t-SNE, UMAP never builds the N² distance matrix — it only keeps these edges.
      </aside>
    </section>

    <!-- 41 · Make the edges fuzzy -->
    <section>
      <div class="eyebrow">Embeddings · UMAP · step 2 of 3</div>
      <h2>Make the edges fuzzy</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('umap_step2_fuzzy.png')" alt="Edges with thickness proportional to their fuzzy weight" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Rescale each edge by local density, so a weight means 'probability this edge exists'.</p>
      <aside class="notes">
        (~30 s) The fuzzy step is the one to dwell on: rescaling by each point's own nearest neighbour makes sparse and dense regions comparable — how UMAP avoids t-SNE's density distortion.
      </aside>
    </section>

    <!-- 42 · Lay the graph out -->
    <section>
      <div class="eyebrow">Embeddings · UMAP · step 3 of 3</div>
      <h2>Lay the graph out</h2>
      <div class="figure" style="display: block; width: 100%; max-width: 780px; margin: 0 auto">
        <img :src="asset('umap_step3_layout.png')" alt="Final UMAP embedding with four separated clusters" style="width: 100%; height: auto" />
      </div>
      <p class="small muted center" style="margin-top: 0.35em; max-width: 760px; margin-inline: auto">Attract along edges, repel sampled non-neighbours — the clusters separate.</p>
      <aside class="notes">
        (~30 s) The layout keeps global structure, unlike t-SNE. But it is still a projector — HDBSCAN* runs afterwards. Bridge to EVoC.
      </aside>
    </section>

    <!-- 43 · EVoC — the fusion -->
    <section>
      <div class="eyebrow">The fusion · EVoC</div>
      <h2>EVoC: embed and cluster in one pass</h2>
      <p class="small">
        <strong>EVōC</strong> — Embedding Vector Oriented Clustering (Tutte Institute; said
        “evoke”) — is not a projector. It is the three ideas we have just built, fused into a
        single fit that hands back <strong>labels, not coordinates</strong>.
      </p>
      <div class="cols" style="--n: 3; margin-top: 0.45em">
        <div class="panel">
          <h3>from UMAP</h3>
          <p class="small">
            The same kNN graph on the raw 16 abundances, laid out by the same node-embedding
            machinery — but into 4–15 dimensions of its own, never a 2-D picture.
          </p>
        </div>
        <div class="panel">
          <h3>from HDBSCAN*</h3>
          <p class="small">
            Mutual-reachability MST and condensed tree on that embedding: density clustering, with
            a noise label for the field stars that belong to nothing.
          </p>
        </div>
        <div class="panel">
          <h3>from PLSCAN</h3>
          <p class="small">
            Persistence over min_cluster_size scores every layer of the hierarchy and returns the
            most persistent one. No scale to guess.
          </p>
        </div>
      </div>
      <div class="cols" style="--n: 1; margin-top: 0.4em">
        <div class="panel flip">
          <p class="small" style="margin: 0">
            <strong>The caveat to hold on to:</strong> EVoC works in <strong>cosine</strong>
            geometry. Our benchmark L2-normalises every abundance vector so Euclidean ≈ cosine —
            that is what keeps the three methods comparable, and it turns out to be the single
            biggest precision lever in the whole study.
          </p>
        </div>
      </div>
      <aside class="notes">
        (~4 min) This is the payoff of the whole tour, so make the three columns explicit: nothing
        in EVoC is new — column one is the UMAP graph, column two is HDBSCAN*, column three is
        PLSCAN. What is new is that they share one fit, so the embedding is built for the
        clustering rather than for your eyes. Correct the obvious misreading before it happens:
        EVoC does not cluster the raw 16-D vectors, it clusters its own internal node embedding —
        the point is that there is no lossy 2-D detour, not that the geometry is untouched. Then
        the caveat: cosine is its native metric, and rather than call that unfair we normalise the
        vectors so all three methods see the same geometry. Next slide opens the hood.
      </aside>
    </section>

    <!-- 44 · EVoC — the pipeline -->
    <section>
      <div class="eyebrow">The fusion · EVoC</div>
      <h2>EVoC, under the hood</h2>
      <ol class="contribs small tight" style="margin-top: 0.2em">
        <li><strong>kNN graph</strong> — n_neighbors neighbours per star, on all 16 abundances, in cosine geometry.</li>
        <li><strong>Node embedding</strong> — UMAP-style layout of that graph, label-propagation init, into 4–15-D.</li>
        <li><strong>Mutual-reachability MST</strong> — Borůvka on the embedding: HDBSCAN*'s own step, unchanged.</li>
        <li><strong>Condensed tree</strong> — the density hierarchy at every scale, exactly as HDBSCAN* builds it.</li>
        <li><strong>Persistence</strong> — a min_cluster_size barcode scores each layer; the most persistent one wins.</li>
      </ol>
      <p class="small muted" style="margin-top: 0.3em">
        Out: <strong>labels + membership strengths + persistence scores + the layer hierarchy</strong> — not one flat partition.
      </p>
      <div class="center" style="margin-top: 0.3em">
        <div class="figure" style="aspect-ratio: 1444 / 443; width: 90%">
          <img :src="asset('evoc_pipeline.png')" alt="EVoC pipeline: kNN graph, node embedding, cluster labels, persistence per layer" style="width: 100%; height: 100%; object-fit: contain" />
        </div>
      </div>
      <aside class="notes">
        (~3 min) Read the five steps as a map back onto the tour: 1–2 is UMAP's graph machinery,
        3–4 is HDBSCAN*, 5 is PLSCAN. Two details are worth pausing on. First, the density
        clustering does not run on the 16-D abundances — it runs on the node embedding, which EVoC
        sizes itself (four dimensions at our n_neighbors=15). Second, step 5 is literally the
        PLSCAN barcode: it sweeps min_cluster_size, scores each resulting layer by total
        persistence, and the winning layer becomes the labels you get back. Point at the four
        panels as you go; the last one is the persistence score per layer.
      </aside>
    </section>

    <!-- 45 · EVoC — why C-space -->
    <section>
      <div class="eyebrow">The fusion · EVoC</div>
      <h2>Why this fits chemical space</h2>
      <div class="slide-body">
      <div class="fig-split" style="--cols: 1.3fr 1fr; margin-top: 0.25em; align-items: start">
        <div>
          <ul class="checklist small">
            <li><strong>Clusters are rare and tiny.</strong> Tens of members hiding in ~10⁵ field stars — you need a method allowed to call most of the data noise. K-means must give every star a home.</li>
            <li><strong>Groups are elongated, and their densities differ wildly.</strong> Compact metal-poor globulars and diffuse open clusters in one catalogue: no single ε or min_cluster_size serves both, so persistence picks the layer instead.</li>
            <li><strong>The signal is the abundance pattern, not its amplitude.</strong> Cosine geometry — EVoC's native metric — is the physically right choice, and it is why we L2-normalise.</li>
            <li><strong>No 2-D detour.</strong> The graph is built on all 16 abundances and the clustering runs in EVoC's own 4–15-D embedding, never in a picture.</li>
          </ul>
          <p class="small muted" style="margin-top: 0.3em">
            Its knobs — noise_level, base_min_cluster_size, n_neighbors, n_epochs — are the ones you
            already met in UMAP and HDBSCAN*.
          </p>
        </div>
        <div class="figure" style="aspect-ratio: 589 / 492">
          <img :src="asset('cspace_corner.png')" alt="Two of the sixteen abundance dimensions: field stars grey, cluster members magenta in one tight clump" style="width: 100%; height: 100%; object-fit: contain" />
        </div>
      </div>
      <p class="small center muted" style="margin-top: 0.35em">
        A plausible story is not evidence. <em>Does the fusion actually recover more clusters on
        APOGEE?</em> That is the benchmark.
      </p>
      </div>
      <aside class="notes">
        (~2 min) Bring the corner plot back from the data slide — it is the shape of the problem in
        two of the sixteen dimensions, and every bullet on the left is a property you can see in it:
        a handful of red points, elongated, sitting inside a grey cloud that is thousands of times
        larger. Each line pairs one property of C-space with the ingredient of EVoC that answers it,
        so this is the slide where the whole narrative arc closes. Then break the spell deliberately:
        this is an argument from assumptions, and assumptions are exactly what the talk has been
        warning about. The benchmark is where the story meets reality.
      </aside>
    </section>

    <!-- 46 · The benchmark -->
    <section>
      <div class="eyebrow">The benchmark · this repo</div>
      <h2>Three methods, one kinematic truth</h2>
      <div class="slide-body">
      <p class="small">
        Reproduce Kos et al. 2017 — cluster by cluster, region by region — then extend it:
        <strong>t-SNE vs UMAP vs EVoC</strong> on APOGEE DR17 + Gaia EDR3, over
        <strong>23 clusters</strong> (16 open, including the Pleiades, and 7 globular;
        Garcia-Dias et al. 2019 + Kos et al. 2017).
      </p>
      <div class="cols" style="--n: 3; margin-top: 0.5em">
        <div class="panel">
          <h3>Pipeline</h3>
          <p class="small">
            allStar DR17 FITS → quality cuts (SNR ≥ 100, ASPCAPFLAG = STARFLAG = 0) → kinematic
            membership labels → 16-D abundance matrix → embed → cluster → score.
          </p>
        </div>
        <div class="panel">
          <h3>Methods</h3>
          <p class="small">
            t-SNE → HDBSCAN* · UMAP → HDBSCAN* · EVoC. The first two only <em>embed</em>, so
            clustering is a second, separate step; EVoC clusters the 16-D vectors directly and
            hands back labels.
          </p>
        </div>
        <div class="panel flip">
          <h3>Score vs the Simbad referee</h3>
          <p class="small">
            <strong>recall</strong> — members recovered by the best-matching predicted cluster.
            <strong>precision</strong> — purity of that cluster. <strong>kNN purity</strong> —
            share of a member's 10 nearest neighbours in the embedding from the same cluster.
            Ground truth is the <strong>Simbad catalogue</strong> (external literature membership),
            kinematic labels as fallback.
          </p>
        </div>
      </div>
      <p class="small muted center" style="margin-top: 0.55em">
        C-space: 15 [X/Fe] + [Fe/H], median-centred and unit-scaled, then row-normalised so Euclidean ≡
        cosine. Stars need ≥ 8 of 16 finite abundances — strict complete-case silently deletes the
        metal-poor globulars.
      </p>
      </div>
      <aside class="notes">
        (~3 min) This is the repo's whole point: the workshop module. t-SNE and UMAP only embed, so
        HDBSCAN* clusters their 2-D maps; EVoC clusters the 16-D vectors internally. Everything is
        scored against the same kinematic ground truth — recall (did we find the members?) and
        precision (were we right?), plus kNN purity as a parameter-free cohesion score. Say out loud
        that kNN purity is the automated version of the paper's "draw a polygon round the group"
        test — no hyperparameters, so no way to tune yourself a good answer. The footnote matters
        too: row-normalisation is what makes Euclidean distance equal cosine, EVoC's native metric,
        and the imputation rule is what keeps M 15 and M 92 in the sample at all.
      </aside>
    </section>

    <!-- 46b · Baseline — re-create the paper -->
    <section>
      <div class="eyebrow">Baseline · re-creating Garcia-Dias et al. 2019</div>
      <h2>First: can we separate the clusters from each other?</h2>
      <p class="small">
        The 2019 paper's question, re-run on our DR17 data: take <strong>only the known cluster
        members</strong> (646 stars, no field), cluster them, and ask how well each star is assigned
        back to its own cluster — scored with the paper's own <strong>homogeneity / v-measure /
        accuracy</strong>.
      </p>
      <div class="cols" style="--n: 2; margin-top: 0.45em">
        <div class="panel">
          <h3>Abundances only</h3>
          <table style="font-size: 0.5em; margin-top: 0.3em">
            <thead>
              <tr><th>method</th><th>homog.</th><th>v-meas.</th><th>acc.</th></tr>
            </thead>
            <tbody>
              <tr><td>t-SNE</td><td>0.26</td><td>0.41</td><td>0.38</td></tr>
              <tr><td>UMAP</td><td>0.32</td><td>0.47</td><td>0.42</td></tr>
              <tr><td>EVoC</td><td><strong>0.50</strong></td><td>0.50</td><td>0.36</td></tr>
            </tbody>
          </table>
          <p class="small muted" style="margin-top: 0.3em">
            Paper reported 0.85 — but with <em>supervised</em> LDA and membership 2σ-clipped in the
            same abundances (circular). Our honest kinematic membership: 0.26–0.50.
          </p>
        </div>
        <div class="panel flip">
          <h3>Abundances + kinematics</h3>
          <table style="font-size: 0.5em; margin-top: 0.3em">
            <thead>
              <tr><th>method</th><th>homog.</th><th>v-meas.</th><th>acc.</th></tr>
            </thead>
            <tbody>
              <tr><td>t-SNE</td><td>0.76</td><td>0.79</td><td>0.72</td></tr>
              <tr><td>UMAP</td><td><strong>0.80</strong></td><td>0.81</td><td>0.75</td></tr>
              <tr><td>EVoC</td><td>0.72</td><td>0.70</td><td>0.55</td></tr>
            </tbody>
          </table>
          <p class="small muted" style="margin-top: 0.3em">
            Add parallax + proper motion + radial velocity: homogeneity jumps to 0.72–0.80. The
            signal chemistry can't find is in the motion.
          </p>
        </div>
      </div>
      <p class="stat center" style="margin-top: 0.45em">Abundances suggest.<br />Kinematics decide.<br />Giants give distance.<br />The main sequence gives age.</p>
      <aside class="notes">
        (~3 min) The baseline the whole talk builds on. Re-run the 2019 question on our data: only
        cluster members, no field, recover which star belongs to which cluster. Two honest
        corrections to the paper's 0.85. First, their best result used LDA, a *supervised*
        projection that already knows the clusters; our numbers are fully unsupervised. Second,
        their membership was 2σ-clipped in the same abundances it then clustered — circular. Our
        kinematic (Gaia) membership has no such leak, and the honest abundances-only number is
        0.26–0.50. Point at t-SNE's completeness 1.0 vs homogeneity 0.26: HDBSCAN merges all the
        similar-age open clusters into one blob — the paper's own "indistinguishable pairs".
        Then the flip panel: add parallax, proper motion and radial velocity and homogeneity
        jumps to 0.72–0.80. Chemistry narrows; kinematics decide. That is the sentence the rest
        of the talk hangs on, and it lands us on the next slide — pulling clusters out of the
        field, not just apart from each other.
      </aside>
    </section>

    <!-- 46c · Confusion matrix -->
    <section>
      <div class="eyebrow">Baseline · the structure behind the score</div>
      <h2>Two families, not twenty clusters</h2>
      <p class="small">
        t-SNE's confusion matrix, row-normalised: each row is a true cluster, each column a
        predicted one. With <strong>abundances alone</strong> every globular collapses into one
        column and every open cluster into another — chemistry sees <em>families</em>, not
        individuals. Add <strong>kinematics</strong> and the diagonal lights up.
      </p>
      <div class="fig-split" style="--cols: 1fr 1.4fr; margin-top: 0.35em; align-items: start">
        <div class="figure">
          <img
            :src="asset('baseline_confusion_chem.png')"
            alt="Confusion matrix for abundances only: two bright columns, one catching all globulars and one all open clusters"
            style="width: 100%; height: auto"
          />
          <div class="caption">
            Abundances only — one globular column (c0), one open-cluster column (c1).
          </div>
        </div>
        <div class="figure">
          <img
            :src="asset('baseline_confusion_kin.png')"
            alt="Confusion matrix for abundances plus kinematics: a bright diagonal recovering the twenty clusters"
            style="width: 100%; height: auto"
          />
          <div class="caption">
            Abundances + kinematics — the diagonal recovers the 20 clusters.
          </div>
        </div>
      </div>
      <aside class="notes">
        (~2 min) The numbers from the last slide, made visible. With abundances alone, t-SNE's
        confusion matrix has essentially two bright columns: c0 catches all six globulars, c1
        catches all fourteen open clusters. That is weak chemical tagging — the paper's own
        conclusion: chemistry identifies the *family* (metal-poor vs solar), not the individual
        cluster. Point at the open-cluster column and name the paper's indistinguishable pairs:
        NGC 2158, NGC 2420 and the Pleiades all share c1. Then the right panel: add parallax,
        proper motion and radial velocity and the diagonal appears — twenty clusters, one column
        each. Same data, same algorithm, one extra feature block. That is the whole story in a
        picture: abundances suggest the family, kinematics decide the cluster.
      </aside>
    </section>

    <!-- 47 · Results -->
    <section>
      <div class="eyebrow">Follow-up · retrieval from the field</div>
      <h2>Chemical tagging is hard — that's the finding</h2>
      <p class="small">
        <strong class="accent">Recall without precision is a blob, not a discovery.</strong>
        UMAP recovers the most members, but its best-matching cluster is ~94% field stars; t-SNE gives
        recall away and buys the only usable precision (0.47 on M 67). Macro kNN purity: 0.11 t-SNE, 0.09 UMAP.
      </p>
      <table style="font-size: 0.52em; margin-top: 0.3em">
        <thead>
          <tr>
            <th>mode</th>
            <th>stars</th>
            <th>t-SNE r / p</th>
            <th>UMAP r / p</th>
            <th>EVoC r / p</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>fast, all-sky</td>
            <td>25k</td>
            <td>0.20 / 0.17</td>
            <td>0.26 / 0.11</td>
            <td>0.54 / 0.01</td>
          </tr>
          <tr>
            <td>region sweep (30°), macro over 22 clusters</td>
            <td>7–27k each</td>
            <td>0.30 / 0.15</td>
            <td>0.55 / 0.06</td>
            <td>0.45 / 0.02</td>
          </tr>
          <tr>
            <td>region, M 67 (30°)</td>
            <td>14k</td>
            <td>0.05 / 0.47</td>
            <td>0.61 / 0.05</td>
            <td>0.49 / 0.05</td>
          </tr>
        </tbody>
      </table>
      <div
        class="figure"
        style="display: block; width: calc(100% - 200px); margin: 0.4em 0 0 200px; padding: 0.35em"
      >
        <img
          :src="asset('benchmark_grid.png')"
          alt="Three panels: t-SNE and UMAP embeddings with true cluster members coloured and field stars grey, and EVoC's labels drawn on the UMAP canvas"
          style="display: block; width: 100%; height: auto"
        />
        <div class="caption">
          Embeddings coloured by true membership, field stars grey. Right panel: EVoC's labels on the
          UMAP canvas — EVoC clusters in 16-D, it never embeds.
        </div>
      </div>
      <aside class="notes">
        (~3 min) Read the table honestly, row by row. Fast all-sky: everything is mediocre. Region
        mode is the paper's own method and it is what moves precision — 0.47 on M 67, the one number
        on this slide you could publish. UMAP's 0.55 recall looks like a win until you read across:
        precision 0.06, so nineteen out of twenty stars in that "cluster" are field. That is the
        blob. Then the picture: same data three ways, members coloured, field grey. Point out that
        the third panel is not an EVoC projection — EVoC has no 2-D output, so its labels are painted
        onto the UMAP canvas. Land it as: chemistry narrows the candidate list, kinematics decide.
      </aside>
    </section>

    <!-- 47b · Isochrone + literature -->
    <section>
      <div class="eyebrow">Follow-up · isochrone fitting vs the literature</div>
      <h2>The same members give back age and distance</h2>
      <p class="small">
        Membership is the input, not the answer. Fit a <strong>PARSEC isochrone</strong> to the
        recovered members (ASteCA + emcee) and read off <strong>age</strong> and
        <strong>distance modulus</strong> — then check against the literature (Dias for open
        clusters, Harris for globulars).
      </p>
      <div class="cols" style="--n: 2; margin-top: 0.45em">
        <div class="panel">
          <h3>Distance modulus — the robust check</h3>
          <table style="font-size: 0.52em; margin-top: 0.3em">
            <thead>
              <tr><th>cluster</th><th>fit dm</th><th>lit dm</th><th>Δ</th></tr>
            </thead>
            <tbody>
              <tr><td>M 67 (APOGEE 2-colour)</td><td>9.65</td><td>9.54</td><td><strong>0.11 ✓</strong></td></tr>
              <tr><td>M 15 (APOGEE 2-colour)</td><td>14.99</td><td>15.04</td><td><strong>0.05 ✓</strong></td></tr>
            </tbody>
          </table>
          <p class="small muted" style="margin-top: 0.3em">
            Two-colour (J−K) photometry pins the distance to ≲ 0.1 mag — a few percent in distance.
          </p>
        </div>
        <div class="panel flip">
          <h3>Age — works, until the grid runs out</h3>
          <table style="font-size: 0.52em; margin-top: 0.3em">
            <thead>
              <tr><th>cluster</th><th>fit age</th><th>lit age</th><th>note</th></tr>
            </thead>
            <tbody>
              <tr><td>M 5 (Gaia CMD)</td><td>11.4</td><td>10.5</td><td>✓</td></tr>
              <tr><td>M 15 (Gaia CMD)</td><td>1.78</td><td>12.5</td><td>✗ grid floor</td></tr>
            </tbody>
          </table>
          <p class="small muted" style="margin-top: 0.3em">
            M 15's [M/H] ≈ −2.4 sits below the metal-poor grid's −2.3 floor; the Gaia-only fallback
            (G − K, G &lt; 16) also leaves the distance degenerate for M 5.
          </p>
        </div>
      </div>
      <p class="stat center" style="margin-top: 0.45em">Distances to 0.1 mag.<br />Ages need the right colour.</p>
      <aside class="notes">
        (~2 min) The bridge from "found the members" to "did the astrophysics". The members from
        the retrieval benchmark feed an isochrone fit: ASteCA's PARSEC grids + an emcee sampler over
        (metallicity, log age, distance modulus, extinction). Two-colour APOGEE photometry (J−K)
        pins the distance modulus to a tenth of a magnitude — M 67 9.65 vs 9.54, M 15 14.99 vs
        15.04 — a few percent in distance, a real astrophysical quantity recovered end-to-end from
        the same abundance-selected stars. Ages are shakier, and that is the honest part to show:
        M 5's Gaia CMD fit gives 11.4 Gyr against 10.5, but M 15 fails — its [M/H] ≈ −2.4 sits
        below the metal-poor grid's −2.3 floor, so the fit walks to the edge and returns a nonsense
        1.8 Gyr. A grid-edge lesson: your result is only as good as the range your models cover.
        If time, point them at the notebook §7–§10 for the HR diagrams, the three-panel
        Kiel/2MASS/Gaia views, and the full literature table.
      </aside>
    </section>

    <!-- 47c · Spectral embeddings -->
    <section>
      <div class="eyebrow">Deep learning · spectral embeddings</div>
      <h2>The RNN latent beats abundances — consistently</h2>
      <p class="small">
        Train a <strong>CNN-LSTM-Attention</strong> network to regress 9 abundances from the raw
        <strong>8575-pixel APOGEE spectrum</strong> (multi-task, RTX 5090, 68k stars across the
        full stellar range — no M-dwarf filter). Cluster on its <strong>256-d latent</strong>, and
        compare to the 16 ASPCAP abundances on the <em>same</em> 878 members and the same scaled
        regions.
      </p>
      <div class="cols" style="--n: 2; margin-top: 0.4em">
        <div class="panel">
          <h3>Cluster-only — homogeneity</h3>
          <table style="font-size: 0.5em; margin-top: 0.25em">
            <thead><tr><th></th><th>t-SNE</th><th>UMAP</th><th>EVoC</th></tr></thead>
            <tbody>
              <tr><td>abundances</td><td>0.25</td><td>0.51</td><td>0.51</td></tr>
              <tr><td><strong>RNN latent (256-d)</strong></td><td><strong>0.40</strong></td><td><strong>0.69</strong></td><td><strong>0.64</strong></td></tr>
            </tbody>
          </table>
        </div>
        <div class="panel flip">
          <h3>Field retrieval — precision (recall)</h3>
          <table style="font-size: 0.5em; margin-top: 0.25em">
            <thead><tr><th></th><th>t-SNE</th><th>UMAP</th><th>EVoC</th></tr></thead>
            <tbody>
              <tr><td>abundances</td><td>0.21 (0.70)</td><td>0.32 (0.54)</td><td>0.28 (0.63)</td></tr>
              <tr><td><strong>RNN latent (256-d)</strong></td><td><strong>0.51</strong> (0.72)</td><td><strong>0.54</strong> (0.68)</td><td><strong>0.56</strong> (0.48)</td></tr>
            </tbody>
          </table>
        </div>
      </div>
      <p class="small muted center" style="margin-top: 0.4em">
        Chemical features only — <strong>no kinematics</strong> (they are the ground truth; adding
        them is circular). The RNN latent is a better chemical-tagging feature than the 16 ASPCAP
        abundances — on both experiments. Globulars tag near-perfectly (M 15 0.97, M 71 0.98,
        M 92 1.00 precision); open clusters stay the hard case (M 67 0.37).
      </p>
      <aside class="notes">
        (~3 min) The deep-learning payoff, now a clean head-to-head. First the setup change that
        matters: the earlier model was trained on M dwarfs only (Teff <= 4100 K), but the cluster
        members are giants — so we rebuilt the training set across the full stellar range (68k
        stars, no filter). Then two experiments, same population, same scaled regions.
        Cluster-only homogeneity: the spectral latent wins every column — t-SNE 0.40 vs 0.25,
        UMAP 0.69 vs 0.51, EVoC 0.64 vs 0.51. Field-retrieval precision: ~2x — 0.51 / 0.54 / 0.56
        vs 0.21 / 0.32 / 0.28 — at comparable recall. Land the two messages. (1) The RNN latent is
        a chemically-meaningful compression of the whole spectrum, not just the 16 ASPCAP elements
        — it sees blends, line wings, and elements ASPCAP skips. (2) Globulars are the showcase
        (M 15, M 71, M 92 near-perfect), open clusters the honest frontier (M 67 0.37) — the same
        "globulars tag cleanly, open clusters don't" lesson, now with a better feature. Close on
        the caveat: the regions are scaled to the cluster (max 3 deg, 10x diameter), so the field
        is modest — a realistic search radius, not a 30-degree blanket.
      </aside>
    </section>

    <!-- 47d · Two surveys -->
    <section>
      <div class="eyebrow">Two surveys &middot; giants + main sequence</div>
      <h2>Distance from the clump, age from the turnoff</h2>
      <p class="small">
        GALAH (DEC &#8818; +25&deg;) sees the <strong>main sequence + turnoff</strong>; APOGEE sees the
        <strong>giants</strong>. Combine them &mdash; 6 clusters, 14 common abundances, Gaia
        parallax/PM/RV &mdash; and the cluster parameters finally separate: the red clump pins the
        distance, the main sequence pins the age.
      </p>
      <div class="cols" style="--n: 2; margin-top: 0.4em">
        <div class="panel">
          <h3>Red-clump distance (J&minus;K)</h3>
          <table style="font-size: 0.5em; margin-top: 0.25em">
            <thead><tr><th>cluster</th><th>dm clump</th><th>dm lit</th><th>resid</th></tr></thead>
            <tbody>
              <tr><td>NGC 6819</td><td>11.92</td><td>11.90</td><td><strong>+0.01</strong></td></tr>
              <tr><td>NGC 2243</td><td>13.08</td><td>13.25</td><td><strong>&minus;0.16</strong></td></tr>
              <tr><td>NGC 7789</td><td>11.68</td><td>11.27</td><td>+0.41</td></tr>
            </tbody>
          </table>
        </div>
        <div class="panel flip">
          <h3>NGC 2243 &mdash; combined fit</h3>
          <table style="font-size: 0.5em; margin-top: 0.25em">
            <thead><tr><th></th><th>age (lit 1.07)</th><th>resid</th><th>std dm</th></tr></thead>
            <tbody>
              <tr><td>APOGEE giants</td><td>1.59</td><td>0.52</td><td>0.111</td></tr>
              <tr><td><strong>MS + giants</strong></td><td><strong>0.98</strong></td><td><strong>0.09</strong></td><td><strong>0.086</strong></td></tr>
            </tbody>
          </table>
        </div>
      </div>
      <p class="small muted center" style="margin-top: 0.4em">
        Chemical tagging is a <strong>two-survey problem</strong>: Kos&rsquo;s GALAH saw the turnoff
        (age), Garcia-Dias&rsquo;s APOGEE saw the giants (distance). Neither alone gives both &mdash;
        together they do.
      </p>
      <aside class="notes">
        (~2 min) The synthesis that closes the loop on the two opening papers. GALAH reaches only
        DEC &lt; +25 deg, so we added two southern intermediate-age clusters (NGC 2243, Collinder 261)
        where both surveys overlap. The red-clump distance is the median K of the J-K clump slice
        (member giants, metallicity-corrected M_K): NGC 6819 lands within 0.01 mag, NGC 2243 within
        0.16 mag. Then the sweet spot: the APOGEE member giants locate the cluster, the GALAH main
        sequence is selected around that kinematic centroid with a relaxed parallax (at 4.7 kpc the
        parallax error rivals the parallax). The combined isochrone fit recovers age 0.98 Gyr
        (resid 0.09 vs 0.52 giants-only) and narrows the distance posterior (std dm 0.111 -> 0.086).
        The honest caveat: the age posterior width is ~tied (std_loga ~0.97) -- age-metallicity
        degeneracy is intrinsic to the CMD, not broken by adding the main sequence. Land: distance
        from giants, age from the main sequence -- one survey for each.
      </aside>
    </section>

    <!-- 48 · Lessons -->
    <section>
      <div class="eyebrow">What the numbers teach</div>
      <h2>Eight lessons from the benchmark</h2>
      <div class="slide-body">
      <ol class="contribs tight small">
        <li><strong>All-sky collapses.</strong> One embedding of every clean star buries each cluster in the field — cut a sky region around it, as Kos et al. did (40° around the Pleiades).</li>
        <li><strong>Precision is the hard part.</strong> Even region by region, macro precision is ≈ 0.15 (t-SNE) and 0.02–0.06 (UMAP, EVoC): the field is chemically similar. Kinematics confirm what chemistry only suggests.</li>
        <li><strong>Recall ≈ 1.0 is usually a blob.</strong> HDBSCAN* swallows the dense field into one cluster, so UMAP "recovers" everything at precision ≈ 0. kNN purity is the honest, parameter-free score.</li>
        <li><strong>Globulars tag cleanly, open clusters don't.</strong> M 3 / M 5 / M 15 reach kNN purity ≈ 0.4–0.5; solar-metallicity open clusters stay ≲ 0.2 — the paper's own caveat.</li>
        <li><strong>t-SNE isolates, UMAP over-merges, EVoC is fast but coarse.</strong> Three sets of assumptions, three answers.</li>
        <li><strong>Element weights are the biggest single lever.</strong> Weighting each element by 1/σ <em>after</em> standardising lifts t-SNE recall 0.09 → 0.44 and EVoC recall 0.38 → 0.45 on M 67 — the noisiest elements were owning every distance.</li>
        <li><strong>Spectra beat abundances.</strong> A CNN-LSTM-Attention latent (256-d) trained on the raw spectrum beats the 16 ASPCAP abundances on every benchmark — and it carries the age signal that element ratios miss (metal-poor globulars: abundances collapse, spectra win).</li>
        <li><strong>Two surveys, two parameters.</strong> APOGEE sees giants (red clump → distance); GALAH sees the main sequence (turnoff → age). Combine them and the isochrone separates — NGC 2243 age residual 0.09.</li>
      </ol>
      <p class="small muted center" style="margin-top: 0.5em">
        And one quieter lesson: row-normalisation is what breaks the blob — precision rises three- to
        fivefold and recall falls. No setting hands you both.
      </p>
      </div>
      <aside class="notes">
        (~3 min) The teaching payload — each lesson maps back to something in the tour. (1) is the
        region detail from the t-SNE slides. (2) is the paper's own "untaggable" caveat, 47 Tuc
        included. (3) is the validation slide's warning about scores that lie, and the reason we
        carry kNN purity at all. (4) is no free lunch in astronomical clothing: the right tool
        depends on the shape of the data, and metal-poor globulars simply are a different shape in
        C-space. (5) is the EVoC caveat about the metric. (6) is the element-weight lever — the concrete knob the hands-on hands them. Close on the footnote: the
        precision/recall trade-off is a knob, not a bug — you choose which error you can live with.
      </aside>
    </section>

    <!-- 49 · Take-away -->
    <section>
      <div class="eyebrow">Take-away</div>
      <h2>No bad models — only mismatched ones</h2>
      <p class="small">
        Every algorithm in this tour encodes a guess: K-means guesses round and equal-sized; DBSCAN
        guesses one density everywhere; HDBSCAN* and PLSCAN drop that guess; t-SNE and UMAP guess that
        the neighbourhood graph is what matters.
      </p>
      <ul class="checklist medium">
        <li>There are no bad models — only models applied outside their assumptions</li>
        <li>Know your data first: shape, scale, density, metric</li>
        <li>Statistics and visualisation matter as much as the model</li>
        <li>Validate against independent evidence — and never quit thinking</li>
      </ul>
      <p class="stat center" style="margin-top: 0.5em">Know the assumptions,<br />know your data</p>
      <aside class="notes">
        (~2 min) Land the through-line one more time, now earned by a concrete astronomy example.
        Walk back up the list on the slide: every clusterer you met today encodes assumptions about
        shape, scale, density and metric, and the failures you saw were never the algorithm being
        "bad" — they were the assumption being wrong for these data. The only way to know you chose
        right is to validate against something the model never saw. Here that was kinematics; in
        their own work it will be something else, but it has to exist.
      </aside>
    </section>

    <!-- 49b · The assignment -->
    <section>
      <div class="eyebrow">The assignment · 25 clusters, ~30 of you</div>
      <h2>Take a cluster, beat our baseline</h2>
      <p class="small">
        The results you just saw are a <strong>baseline, not a ceiling</strong>. Each of you gets one
        cluster from the region sweep — <code>docs/region_sweep_results.md</code> in the repo holds our
        current recall / precision / kNN-purity for every one. Your job: <strong>extend the pipeline and
        get better numbers for your cluster</strong>.
      </p>
      <div class="cols" style="--n: 3; margin-top: 0.45em">
        <div class="panel">
          <h3>Your cluster</h3>
          <p class="small">
            18 open clusters (Pleiades, M 67, NGC 188, NGC 6819, NGC 2243, Collinder 261, …) and 7 globulars (M 3, M 5, M 13,
            M 15, M 71, M 107, M 92). Pair up on the crowded fields. The globulars already tag
            cleanly — your real job is the open clusters that don't.
          </p>
        </div>
        <div class="panel">
          <h3>Levers we've already found</h3>
          <p class="small">
            Element-wise 1/σ weights (after standardising) — the single biggest gain. HDBSCAN
            <code>min_cluster_size</code>, t-SNE <code>perplexity</code>, SNR cut, dwarf-only,
            element subset. Newer levers: the <strong>RNN spectral latent</strong> (256-d), the
            <strong>red-clump / isochrone</strong> fit, and the <strong>GALAH main-sequence</strong> cross-match.
            All logged in <code>docs/experiment_results.md</code>.
          </p>
        </div>
        <div class="panel flip">
          <h3>Deliverable</h3>
          <p class="small">
            One slide: your cluster, your best recall / precision / purity, and the <em>one</em> change
            that moved them. Beat the baseline — or explain honestly why your cluster resists. Either
            answer is a result.
          </p>
        </div>
      </div>
      <p class="stat center" style="margin-top: 0.5em">Not a bigger number —<br />knowing which knob moved why</p>
      <aside class="notes">
        (~2 min) This is the hinge between "here is what I did" and "now it is your turn". Hand out
        the cluster list, then make the framing explicit: the table in
        docs/region_sweep_results.md is the scoreboard, and the levers in
        docs/experiment_results.md are the map — we already found the biggest one (element weights,
        0.09 → 0.44 t-SNE recall on M 67), so the easy gains are spoken for. Their real job is the
        next lever, specific to their cluster. Stress that a null result is still a result: if a
        cluster refuses to tag cleanly, say why — that is the paper's own finding for 47 Tuc. Close
        by pointing at the next slide, the three commands that get them running.
      </aside>
    </section>

    <!-- 50 · Hands-on -->
    <section class="title-slide center">
      <div class="eyebrow">Hands-on module</div>
      <h1>Run it yourself</h1>
      <p class="subtitle">
        Reproduce and extend Kos et al. 2017 — t-SNE vs UMAP vs EVoC — in ~1 minute, on a laptop or in Docker.
      </p>
      <pre class="small" style="text-align: left; max-width: 28em; margin: 0.6em auto"><code>uv sync &amp;&amp; uv run cluster download
uv run cluster run --fast
uv run marimo edit notebooks/chemical_tagging.py</code></pre>
      <p class="byline">
        <a href="https://github.com/garciadias/iaa-advanced-neural-networks-2026-draft" target="_blank" rel="noopener">github.com/garciadias/iaa-advanced-neural-networks-2026-draft</a>
      </p>
      <p class="small muted" style="margin-top: 0.35em">
        Full walkthrough &middot; <code>docs/student_activities.md</code> &nbsp;&middot;&nbsp;
        cluster assignment &middot; <code>docs/cluster_assignment.md</code> &nbsp;&middot;&nbsp;
        Docker-only path &middot; <code>docs/docker.md</code>
      </p>
      <div class="cols" style="--n: 3; gap: 0.7em; max-width: 22em; margin: 0.8em auto 0; align-items: stretch">
        <div>
          <p class="figure-placeholder__desc" style="margin-bottom: 0.35em">Workshop repo</p>
          <div class="figure" style="display: block; aspect-ratio: 1 / 1; padding: 0.4em">
            <img
              :src="asset('presentation.png')"
              alt="QR code linking to the workshop repository and these slides"
              title="Workshop repository"
              style="display: block; width: 100%; height: 100%; object-fit: contain"
            />
          </div>
        </div>
        <div>
          <p class="figure-placeholder__desc" style="margin-bottom: 0.35em">Flip a flag</p>
          <div class="panel" style="aspect-ratio: 1 / 1; min-height: 0; padding: 0.55em; display: flex; align-items: center; justify-content: center">
            <code style="font-size: 0.5em; text-align: left; line-height: 1.5">--fast<br />--full<br />--cluster "M 67"<br />--region 30</code>
          </div>
        </div>
        <div>
          <p class="figure-placeholder__desc" style="margin-bottom: 0.35em">Three methods</p>
          <div class="panel flip" style="aspect-ratio: 1 / 1; min-height: 0; padding: 0.55em; display: flex; align-items: center; justify-content: center">
            <span style="font-size: 0.75em; text-align: center">t-SNE · UMAP · EVoC</span>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~2 min) Close on the workshop. Three commands: install and pull the allStar file (3.7 GB, so
        do it on the hotel wifi tonight, not now), one run, one notebook. The fast run caps the field
        at 25 000 stars and finishes in about a minute; `--full` drops the cap and takes ten to
        twenty. Invite them to add `--cluster "M 67" --region 30` and watch the precision column jump
        — that is lesson one from the previous slide, reproduced on their own laptop in sixty seconds.
      </aside>
    </section>

    <!-- 51 · References -->
    <section>
      <div class="eyebrow">References</div>
      <h2>Papers &amp; code behind this talk</h2>
      <div class="cols" style="--n: 3; margin-top: 0.3em; align-items: start">
        <div class="panel">
          <h3>The science &amp; this work</h3>
          <ul style="list-style: none; margin: 0; font-size: 0.5em; line-height: 1.1">
            <li style="margin-bottom: 0.12em"><a href="https://arxiv.org/abs/1709.00794" target="_blank" rel="noopener">Kos et al. 2017</a> — GALAH: chemical tagging of star clusters &amp; new members in the Pleiades. <span class="muted">arXiv:1709.00794</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://ui.adsabs.harvard.edu/abs/2002ARA%26A..40..487F/abstract" target="_blank" rel="noopener">Freeman &amp; Bland-Hawthorn 2002</a> — The New Galaxy: signatures of its formation. <span class="muted">ARA&amp;A 40, 487</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://doi.org/10.1051/0004-6361/201935223" target="_blank" rel="noopener">Garcia-Dias et al. 2019</a> — Machine learning in APOGEE: stellar populations. <span class="muted">A&amp;A 629, A34</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://doi.org/10.1016/B978-0-12-815739-8.00013-4" target="_blank" rel="noopener">Garcia-Dias et al. 2020</a> — Clustering analysis (book chapter). <span class="muted">Elsevier</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://arxiv.org/abs/2512.16558" target="_blank" rel="noopener">Bot, McInnes &amp; Aerts 2025</a> — Persistent multiscale density-based clustering (PLSCAN). <span class="muted">arXiv:2512.16558</span></li>
          </ul>
        </div>
        <div class="panel">
          <h3>Clustering &amp; density algorithms</h3>
          <ul style="list-style: none; margin: 0; font-size: 0.5em; line-height: 1.1">
            <li style="margin-bottom: 0.12em"><a href="https://doi.org/10.1109/TIT.1982.1056489" target="_blank" rel="noopener">Lloyd 1982</a> — Least squares quantization in PCM (K-means). <span class="muted">IEEE TIT</span></li>
            <li style="margin-bottom: 0.12em">MacQueen 1967 — Some methods for classification &amp; analysis of multivariate observations (K-means). <span class="muted">Berkeley Symp.</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf" target="_blank" rel="noopener">Arthur &amp; Vassilvitskii 2007</a> — k-means++: the advantages of careful seeding. <span class="muted">SODA</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://www.aaai.org/Papers/KDD/1996/KDD96-037.pdf" target="_blank" rel="noopener">Ester et al. 1996</a> — A density-based algorithm for discovering clusters (DBSCAN). <span class="muted">KDD</span></li>
            <li style="margin-bottom: 0.12em">Ankerst et al. 1999 — OPTICS: ordering points to identify the clustering structure. <span class="muted">SIGMOD</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://doi.org/10.1007/978-3-642-37456-2_14" target="_blank" rel="noopener">Campello et al. 2013</a> — Density-based clustering via hierarchical density estimates. <span class="muted">PAKDD</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://doi.org/10.1145/2733381" target="_blank" rel="noopener">Campello et al. 2015</a> — Hierarchical density estimates (HDBSCAN*). <span class="muted">ACM TKDD</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://arxiv.org/abs/1705.07321" target="_blank" rel="noopener">McInnes &amp; Healy 2017</a> — Accelerated hierarchical density clustering (the <code>hdbscan</code> library). <span class="muted">arXiv:1705.07321</span></li>
          </ul>
        </div>
        <div class="panel">
          <h3>Embeddings &amp; validation</h3>
          <ul style="list-style: none; margin: 0; font-size: 0.5em; line-height: 1.1">
            <li style="margin-bottom: 0.12em"><a href="https://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf" target="_blank" rel="noopener">van der Maaten &amp; Hinton 2008</a> — Visualizing data using t-SNE. <span class="muted">JMLR 9</span></li>
            <li style="margin-bottom: 0.12em">van der Maaten 2014 — Accelerating t-SNE using tree-based algorithms (Barnes-Hut). <span class="muted">JMLR 15</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://arxiv.org/abs/1802.03426" target="_blank" rel="noopener">McInnes, Healy &amp; Melville 2018</a> — UMAP: uniform manifold approximation &amp; projection. <span class="muted">arXiv:1802.03426</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://doi.org/10.1016/0377-0427(87)90125-7" target="_blank" rel="noopener">Rousseeuw 1987</a> — Silhouettes: a graphical aid. <span class="muted">J. Comput. Appl. Math.</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://doi.org/10.1214/aos/1176346577" target="_blank" rel="noopener">Hartigan &amp; Hartigan 1985</a> — The dip test of unimodality. <span class="muted">Ann. Stat.</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://aclanthology.org/D07-1043.pdf" target="_blank" rel="noopener">Rosenberg &amp; Hirschberg 2007</a> — V-measure: external cluster evaluation. <span class="muted">EMNLP</span></li>
            <li style="margin-bottom: 0.12em"><a href="https://doi.org/10.1016/j.neuroimage.2009.06.014" target="_blank" rel="noopener">Nanetti et al. 2009</a> — Repeated K-means cortical parcellation. <span class="muted">NeuroImage</span></li>
          </ul>
        </div>
      </div>
      <p class="small" style="margin-top: 0.45em">
        Code:
        <a href="https://github.com/TutteInstitute/evoc" target="_blank" rel="noopener">EVoC</a> ·
        <a href="https://github.com/jelmerbot/fast_plscan" target="_blank" rel="noopener">PLSCAN</a> ·
        <a href="https://github.com/garciadias/iaa-advanced-neural-networks-2026-draft" target="_blank" rel="noopener">workshop repo</a>
      </p>
      <aside class="notes">
        (~1 min)
        Appendix — do not present line-by-line. Leave it up while taking questions, or point people
        here for the arXiv/DOI links. Every claim in the talk traces to one of these. Two entries
        carry no link on purpose: MacQueen 1967 (proceedings, no stable DOI) and van der Maaten 2014,
        where the JMLR volume needs checking before a URL goes on the slide.
      </aside>
    </section>
  </RevealDeck>
</template>

