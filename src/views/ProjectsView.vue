<script setup>
import projects from '@/content/projects.json'
</script>

<template>
  <div class="container-wide py-12 sm:py-16">
    <div class="section-eyebrow">Portfolio</div>
    <h1 class="text-3xl sm:text-4xl font-bold text-gray-900 dark:text-white mb-10">Projects</h1>

    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <article
        v-for="p in projects"
        :key="p.id"
        class="card p-6 flex flex-col"
      >
        <div class="flex items-baseline justify-between gap-2">
          <h2 class="font-display font-bold text-xl text-gray-900 dark:text-white">{{ p.title }}</h2>
          <span class="text-xs font-mono text-gray-500 dark:text-gray-400 whitespace-nowrap">{{ p.year }}</span>
        </div>
        <p class="mt-1 text-sm font-medium text-primary-600 dark:text-primary-400">{{ p.tagline }}</p>
        <p class="mt-3 text-sm text-gray-700 dark:text-gray-300 leading-relaxed flex-1">{{ p.description }}</p>

        <div class="mt-4 flex flex-wrap gap-1.5">
          <span v-for="t in p.tags" :key="t" class="chip">{{ t }}</span>
        </div>

        <div class="mt-5 flex flex-wrap gap-3">
          <a
            v-for="l in p.links"
            :key="l.href"
            :href="l.href"
            target="_blank"
            rel="noopener"
            class="text-sm font-semibold link inline-flex items-center gap-1"
          >
            {{ l.label }}
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" class="w-3.5 h-3.5"><path d="M7 17L17 7M9 7h8v8" stroke-linecap="round" stroke-linejoin="round"/></svg>
          </a>
        </div>
      </article>
    </div>
  </div>
</template>
