<script setup>
// Secret full index at /presentations/list — shows ALL decks, including those
// flagged `unlisted` that are hidden from the public /presentations gallery.
// This route is intentionally not linked from anywhere in the site nav.
import { presentations } from '@/content/presentations'
</script>

<template>
  <div class="container-wide py-12 sm:py-16">
    <div class="section-eyebrow">Full index</div>
    <h1 class="text-3xl sm:text-4xl font-bold text-gray-900 dark:text-white mb-2">
      All presentations
    </h1>
    <p class="text-gray-600 dark:text-gray-400 mb-10 max-w-2xl">
      Every deck, including unlisted ones not shown on the public
      <router-link
        :to="{ name: 'presentations' }"
        class="underline decoration-dotted hover:text-primary-500"
        >presentations page</router-link
      >. Share these links directly.
    </p>

    <div
      v-if="presentations.length === 0"
      class="card p-10 text-center text-gray-500 dark:text-gray-400"
    >
      No presentations published yet — check back soon.
    </div>

    <div v-else class="grid grid-cols-1 md:grid-cols-2 gap-8">
      <router-link
        v-for="p in presentations"
        :key="p.id"
        :to="{ name: 'presentation', params: { id: p.id } }"
        class="card group p-0 overflow-hidden flex flex-col"
      >
        <!-- Deck preview banner — slide cover image when present, else gradient -->
        <div
          class="aspect-video relative flex flex-col justify-end p-6 bg-gradient-to-br from-primary-600 via-primary-700 to-dracula-bg overflow-hidden"
        >
          <template v-if="p.cover">
            <img
              :src="p.cover"
              :alt="`${p.title} — cover`"
              class="absolute inset-0 w-full h-full object-cover"
            />
            <!-- scrim so the title stays legible over the image -->
            <div
              class="absolute inset-0 bg-gradient-to-t from-black/90 via-black/55 to-black/15"
            ></div>
          </template>
          <div v-else class="absolute inset-0 bg-grid opacity-30"></div>
          <div class="relative">
            <div
              class="text-[10px] font-mono uppercase tracking-widest text-white/70 mb-1"
            >
              {{ p.venue }} · {{ p.date }}
            </div>
            <h2 class="font-display font-bold text-xl text-white leading-tight drop-shadow">
              {{ p.title }}
            </h2>
            <p v-if="p.subtitle" class="text-sm text-white/80 mt-1 drop-shadow">{{ p.subtitle }}</p>
          </div>

          <!-- Unlisted marker -->
          <span
            v-if="p.unlisted"
            class="absolute top-4 left-4 inline-flex items-center gap-1.5 rounded-full bg-amber-500/85 backdrop-blur px-3 py-1 text-xs font-medium text-black ring-1 ring-amber-200/40"
          >
            <svg
              xmlns="http://www.w3.org/2000/svg"
              viewBox="0 0 24 24"
              fill="none"
              stroke="currentColor"
              stroke-width="2"
              stroke-linecap="round"
              stroke-linejoin="round"
              class="w-3.5 h-3.5"
            >
              <path d="M17.94 17.94A10.07 10.07 0 0 1 12 20c-7 0-11-8-11-8a18.45 18.45 0 0 1 5.06-5.94M9.9 4.24A9.12 9.12 0 0 1 12 4c7 0 11 8 11 8a18.5 18.5 0 0 1-2.16 3.19m-6.72-1.07a3 3 0 1 1-4.24-4.24" />
              <line x1="1" y1="1" x2="23" y2="23" />
            </svg>
            Unlisted
          </span>

          <span
            class="absolute top-4 right-4 inline-flex items-center gap-1.5 rounded-full bg-white/15 backdrop-blur px-3 py-1 text-xs font-medium text-white ring-1 ring-white/25 group-hover:bg-white/25 transition"
          >
            <svg
              xmlns="http://www.w3.org/2000/svg"
              viewBox="0 0 24 24"
              fill="currentColor"
              class="w-3.5 h-3.5"
            >
              <path d="M8 5v14l11-7z" />
            </svg>
            Open deck
          </span>
        </div>

        <div class="p-6 flex flex-col flex-1">
          <p class="text-sm text-gray-700 dark:text-gray-300 leading-relaxed flex-1">
            {{ p.description }}
          </p>
          <div class="mt-4 flex flex-wrap gap-1.5">
            <span v-for="tag in p.tags" :key="tag" class="chip">{{ tag }}</span>
          </div>
        </div>
      </router-link>
    </div>
  </div>
</template>
