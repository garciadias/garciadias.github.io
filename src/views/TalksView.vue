<script setup>
import talks from '@/content/talks.json'
</script>

<template>
  <div class="container-wide py-12 sm:py-16">
    <div class="section-eyebrow">Speaking</div>
    <h1 class="text-3xl sm:text-4xl font-bold text-gray-900 dark:text-white mb-2">Talks & Videos</h1>
    <p class="text-gray-600 dark:text-gray-400 mb-10 max-w-2xl">
      Recordings of talks, demos and conference appearances. New content is added by editing
      <code class="font-mono text-xs px-1.5 py-0.5 rounded bg-gray-100 dark:bg-gray-800">src/content/talks.json</code>.
    </p>

    <div v-if="talks.length === 0" class="card p-10 text-center text-gray-500 dark:text-gray-400">
      No talks published yet — check back soon.
    </div>

    <div v-else class="grid grid-cols-1 md:grid-cols-2 gap-8">
      <article v-for="t in talks" :key="t.id" class="card overflow-hidden flex flex-col">
        <div class="aspect-video bg-gray-100 dark:bg-gray-800">
          <iframe
            v-if="t.youtubeId"
            :src="`https://www.youtube-nocookie.com/embed/${t.youtubeId}`"
            :title="t.title"
            loading="lazy"
            allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
            allowfullscreen
            class="w-full h-full"
          ></iframe>
        </div>
        <div class="p-6 flex flex-col flex-1">
          <div class="flex items-baseline justify-between gap-2">
            <h2 class="font-display font-bold text-lg text-gray-900 dark:text-white">{{ t.title }}</h2>
            <span class="text-xs font-mono text-gray-500 dark:text-gray-400">{{ t.year }}</span>
          </div>
          <div class="mt-1 text-sm text-gray-600 dark:text-gray-400">
            {{ t.speaker }} · {{ t.venue }}
          </div>
          <p class="mt-3 text-sm text-gray-700 dark:text-gray-300 leading-relaxed flex-1">
            {{ t.description }}
          </p>
          <div class="mt-4 flex flex-wrap gap-1.5">
            <span v-for="tag in t.tags" :key="tag" class="chip">{{ tag }}</span>
          </div>
        </div>
      </article>
    </div>
  </div>
</template>
