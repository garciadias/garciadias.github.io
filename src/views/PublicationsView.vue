<script setup>
import { computed } from 'vue'
import publications from '@/content/publications.json'
import { profile } from '@/content/profile'

const sorted = computed(() => [...publications].sort((a, b) => b.year - a.year))
</script>

<template>
  <div class="container-prose py-12 sm:py-16">
    <div class="section-eyebrow">Research output</div>
    <h1 class="text-3xl sm:text-4xl font-bold text-gray-900 dark:text-white mb-2">Publications</h1>
    <p class="text-gray-600 dark:text-gray-400 mb-10">
      18 Q1 journal articles · h-index 19 · 5,000+ citations.
      Full list on
      <a
        href="https://scholar.google.co.uk/citations?hl=en&user=MlwZerQAAAAJ"
        target="_blank"
        rel="noopener"
        class="link"
        >Google Scholar</a
      >. Selected highlights below.
    </p>

    <ul class="space-y-4">
      <li v-for="(p, i) in sorted" :key="i" class="card p-6">
        <div class="flex flex-wrap items-baseline justify-between gap-2">
          <h2 class="font-display font-semibold text-lg text-gray-900 dark:text-white leading-snug">
            {{ p.title }}
          </h2>
          <span class="text-xs font-mono text-gray-500 dark:text-gray-400">{{ p.year }}</span>
        </div>
        <div class="mt-1 text-sm text-gray-600 dark:text-gray-400">
          <span class="font-medium">{{ p.venue }}</span>
          <span class="mx-1.5">·</span>
          <span>{{ p.type }}</span>
        </div>
        <p v-if="p.description" class="mt-3 text-sm text-gray-700 dark:text-gray-300 leading-relaxed">
          {{ p.description }}
        </p>
        <div class="mt-3 flex flex-wrap gap-3">
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
      </li>
    </ul>
  </div>
</template>
