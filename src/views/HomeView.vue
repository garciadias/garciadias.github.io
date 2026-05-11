<script setup>
import { profile } from '@/content/profile'
</script>

<template>
  <!-- Hero -->
  <section class="relative overflow-hidden">
    <div class="absolute inset-0 bg-grid pointer-events-none" aria-hidden="true"></div>
    <div
      class="absolute inset-x-0 top-0 -z-10 h-[400px] bg-gradient-to-b from-primary-50 via-white to-white dark:from-primary-950/40 dark:via-gray-950 dark:to-gray-950"
      aria-hidden="true"
    ></div>

    <div class="container-wide relative pt-16 pb-20 sm:pt-24 sm:pb-28">
      <div class="grid grid-cols-1 lg:grid-cols-12 gap-10 items-center">
        <div class="lg:col-span-8">
          <div class="inline-flex items-center gap-2 chip mb-6">
            <span class="w-2 h-2 rounded-full bg-emerald-500 animate-pulse"></span>
            Available for collaboration · {{ profile.location }}
          </div>
          <h1 class="text-4xl sm:text-5xl lg:text-6xl font-extrabold tracking-tight text-gray-900 dark:text-white">
            {{ profile.name }}
          </h1>
          <p class="mt-4 text-lg sm:text-xl text-gray-600 dark:text-gray-300 max-w-3xl leading-relaxed">
            {{ profile.headline }}
          </p>
          <p class="mt-6 text-base sm:text-lg text-gray-700 dark:text-gray-400 max-w-3xl leading-relaxed">
            {{ profile.intro }}
          </p>
          <div class="mt-8 flex flex-wrap gap-3">
            <router-link
              :to="{ name: 'projects' }"
              class="inline-flex items-center gap-2 rounded-lg bg-primary-600 hover:bg-primary-700 text-white px-5 py-2.5 text-sm font-semibold shadow-sm transition"
            >
              See projects
              <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" class="w-4 h-4"><path d="M5 12h14M13 6l6 6-6 6" stroke-linecap="round" stroke-linejoin="round"/></svg>
            </router-link>
          </div>
          <p class="mt-4 text-xs text-gray-500 dark:text-gray-500">{{ profile.rightToWork }}</p>
        </div>

        <div class="lg:col-span-4">
          <div class="grid grid-cols-2 gap-3">
            <img src="/static/img/3.png" alt="Rafael with the Neurofind team" class="rounded-xl object-cover w-full h-32 sm:h-40 col-span-2" />
            <img src="/static/img/5.png" alt="Portrait" class="rounded-xl object-cover w-full h-32 sm:h-40" />
            <img src="/static/img/4.png" alt="Giving a talk about GitHub" class="rounded-xl object-cover w-full h-32 sm:h-40" />
          </div>
        </div>
      </div>
    </div>
  </section>

  <!-- Highlights -->
  <section class="container-wide py-12">
    <div class="section-eyebrow">Selected impact</div>
    <h2 class="section-title mb-8">Things I’ve helped ship</h2>
    <div class="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
      <a
        v-for="h in profile.highlights"
        :key="h.title"
        :href="h.href"
        target="_blank"
        rel="noopener"
        class="card p-6 group"
      >
        <div class="flex items-start justify-between gap-3">
          <h3 class="font-display font-bold text-lg text-gray-900 dark:text-white group-hover:text-primary-600 dark:group-hover:text-primary-400">
            {{ h.title }}
          </h3>
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" class="w-4 h-4 text-gray-400 group-hover:text-primary-500 mt-1.5"><path d="M7 17L17 7M9 7h8v8" stroke-linecap="round" stroke-linejoin="round"/></svg>
        </div>
        <p class="mt-2 text-sm text-gray-600 dark:text-gray-400">{{ h.description }}</p>
      </a>
    </div>
  </section>

  <!-- Stack -->
  <section class="container-wide py-12">
    <div class="section-eyebrow">What I work with</div>
    <h2 class="section-title mb-8">Core stack</h2>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div
        v-for="(items, group) in profile.stack"
        :key="group"
        class="card p-6"
      >
        <h3 class="font-mono text-xs uppercase tracking-widest text-primary-600 dark:text-primary-400 mb-3">
          {{ group }}
        </h3>
        <div class="flex flex-wrap gap-2">
          <span v-for="item in items" :key="item" class="chip">{{ item }}</span>
        </div>
      </div>
    </div>
  </section>
</template>
