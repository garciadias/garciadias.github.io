<script setup>
import { experience, education } from '@/content/experience'
import { profile } from '@/content/profile'
</script>

<template>
  <div class="container-prose py-12 sm:py-16">
    <div class="section-eyebrow">Experience</div>
    <div class="flex items-end justify-between flex-wrap gap-4 mb-10">
      <h1 class="text-3xl sm:text-4xl font-bold text-gray-900 dark:text-white">CV & Experience</h1>
      <a
        :href="profile.cvPath"
        target="_blank"
        rel="noopener"
        class="inline-flex items-center gap-2 rounded-lg border border-gray-300 dark:border-gray-700 px-4 py-2 text-sm font-semibold hover:bg-gray-50 dark:hover:bg-gray-800"
      >
        Download PDF
      </a>
    </div>

    <ol class="relative border-l border-gray-200 dark:border-gray-800 space-y-10 pl-6">
      <li v-for="(role, idx) in experience" :key="idx" class="relative">
        <span class="absolute -left-[31px] top-1.5 w-3 h-3 rounded-full bg-primary-500 ring-4 ring-white dark:ring-gray-950"></span>

        <div class="card p-6">
          <div class="flex flex-wrap items-baseline justify-between gap-2">
            <h2 class="font-display font-bold text-xl text-gray-900 dark:text-white">{{ role.role }}</h2>
            <span class="text-xs font-mono text-gray-500 dark:text-gray-400">{{ role.period }}</span>
          </div>
          <div class="mt-1 text-sm text-gray-600 dark:text-gray-400">
            {{ role.org }} · {{ role.location }}
          </div>
          <p v-if="role.summary" class="mt-3 text-sm text-gray-700 dark:text-gray-300 leading-relaxed">
            {{ role.summary }}
          </p>

          <div v-for="(s, i) in role.sections" :key="i" class="mt-5">
            <h3 v-if="s.heading" class="font-semibold text-gray-900 dark:text-white">
              <a v-if="s.href" :href="s.href" target="_blank" rel="noopener" class="link">{{ s.heading }}</a>
              <span v-else>{{ s.heading }}</span>
              <span v-if="s.subheading" class="ml-2 text-xs font-normal text-gray-500 dark:text-gray-400">{{ s.subheading }}</span>
            </h3>
            <ul class="mt-2 space-y-2 text-sm text-gray-700 dark:text-gray-300 leading-relaxed list-disc pl-5 marker:text-primary-500">
              <li v-for="(b, bi) in s.bullets" :key="bi">{{ b }}</li>
            </ul>
          </div>
        </div>
      </li>
    </ol>

    <section class="mt-16">
      <div class="section-eyebrow">Education</div>
      <h2 class="section-title mb-6">Academic background</h2>
      <div class="space-y-3">
        <div v-for="(e, i) in education" :key="i" class="card p-5 flex flex-wrap items-baseline justify-between gap-2">
          <div>
            <div class="font-semibold text-gray-900 dark:text-white">{{ e.degree }}</div>
            <div class="text-sm text-gray-600 dark:text-gray-400">{{ e.institution }}</div>
          </div>
          <span class="text-xs font-mono text-gray-500 dark:text-gray-400">{{ e.year }}</span>
        </div>
      </div>
    </section>

    <section class="mt-12">
      <div class="section-eyebrow">Languages</div>
      <div class="flex flex-wrap gap-2">
        <span v-for="l in profile.languages" :key="l" class="chip">{{ l }}</span>
      </div>
    </section>
  </div>
</template>
