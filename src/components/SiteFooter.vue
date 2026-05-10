<script setup>
import { profile } from '@/content/profile'
import SocialIcon from '@/components/SocialIcon.vue'

const year = new Date().getFullYear()
</script>

<template>
  <footer class="border-t border-gray-200 dark:border-gray-800 mt-16">
    <div class="container-wide py-8 flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4">
      <div class="text-sm text-gray-500 dark:text-gray-400">
        © {{ year }} {{ profile.shortName }} ·
        <a
          href="http://creativecommons.org/licenses/by-sa/4.0/"
          class="link"
          target="_blank"
          rel="noopener"
          >CC BY-SA 4.0</a
        >
      </div>
      <div class="flex items-center gap-1">
        <a
          v-for="s in profile.socials.filter((s) => s.icon !== 'mail')"
          :key="s.label"
          :href="s.href"
          :aria-label="s.label"
          target="_blank"
          rel="noopener"
          class="p-2 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-900 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-800"
        >
          <SocialIcon :name="s.icon" />
        </a>
      </div>
    </div>
  </footer>
</template>
