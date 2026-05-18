<script setup>
import { experience, education } from '@/content/experience'
import { profile } from '@/content/profile'
import { generateCVPDF } from '@/utils/generateCVPDF'

const downloadCV = () => {
  generateCVPDF(experience, education, profile)
}
</script>

<template>
  <div class="container-prose py-12 sm:py-16">
    <div class="section-eyebrow">Working life</div>
    <div class="flex items-end justify-between flex-wrap gap-4 mb-4">
      <h1 class="text-3xl sm:text-4xl font-bold text-gray-900 dark:text-white">Career</h1>
      <button
        @click="downloadCV"
        class="inline-flex items-center gap-2 rounded-lg border border-gray-300 dark:border-gray-700 px-4 py-2 text-sm font-semibold hover:bg-gray-50 dark:hover:bg-gray-800 cursor-pointer"
      >
        Download CV (PDF)
      </button>
    </div>
    <p class="text-gray-600 dark:text-gray-400 mb-10 max-w-2xl">
      A decade of work spanning healthcare AI, foundation models, credit risk and astrophysics.
    </p>

    <div v-if="experience && experience.length" class="space-y-6">
      <article
        v-for="(role, idx) in experience"
        :key="idx"
        class="card p-6"
      >
        <div class="flex flex-wrap items-baseline justify-between gap-2">
          <h2 class="font-display font-bold text-xl text-gray-900 dark:text-white">{{ role.role }}</h2>
          <span class="text-xs font-mono text-gray-500 dark:text-gray-400 whitespace-nowrap">{{ role.period }}</span>
        </div>
        <div class="mt-1 text-sm text-primary-600 dark:text-primary-400 font-medium">
          {{ role.org }}<span v-if="role.location"> · {{ role.location }}</span>
        </div>
        <p v-if="role.summary" class="mt-3 text-sm text-gray-700 dark:text-gray-300 leading-relaxed">
          {{ role.summary }}
        </p>

        <template v-for="(s, i) in role.sections" :key="i">
          <div class="mt-5">
            <h3 v-if="s.heading" class="font-semibold text-gray-900 dark:text-white">
              <a v-if="s.href" :href="s.href" target="_blank" rel="noopener" class="link">{{ s.heading }}</a>
              <span v-else>{{ s.heading }}</span>
              <span v-if="s.subheading" class="ml-2 text-xs font-normal text-gray-500 dark:text-gray-400">{{ s.subheading }}</span>
            </h3>
            <ul v-if="s.bullets && s.bullets.length" class="mt-2 space-y-2 text-sm text-gray-700 dark:text-gray-300 leading-relaxed list-disc pl-5">
              <li v-for="(b, bi) in s.bullets" :key="bi">{{ b }}</li>
            </ul>
          </div>
        </template>
      </article>
    </div>
    <p v-else class="text-gray-500 dark:text-gray-400">No experience entries available.</p>

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
