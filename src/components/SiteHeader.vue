<script setup>
import { ref } from 'vue'
import { useDarkMode } from '@/composables/useDarkMode'
import { profile } from '@/content/profile'
import SocialIcon from '@/components/SocialIcon.vue'

const { isDark, toggle } = useDarkMode()
const mobileOpen = ref(false)

const links = [
  { name: 'home', label: 'Home' },
  { name: 'cv', label: 'Experience' },
  { name: 'projects', label: 'Projects' },
  { name: 'publications', label: 'Publications' },
  { name: 'talks', label: 'Talks' }
]
</script>

<template>
  <header
    class="sticky top-0 z-30 backdrop-blur bg-white/80 dark:bg-gray-950/80 border-b border-gray-200/80 dark:border-gray-800/80"
  >
    <div class="container-wide flex items-center justify-between h-16">
      <router-link
        :to="{ name: 'home' }"
        class="font-display font-bold text-lg sm:text-xl tracking-tight text-gray-900 dark:text-white hover:text-primary-600 dark:hover:text-primary-400"
      >
        Rafael <span class="text-primary-600 dark:text-primary-400">Garcia-Dias</span>
      </router-link>

      <nav class="hidden md:flex items-center gap-1">
        <router-link
          v-for="l in links"
          :key="l.name"
          :to="{ name: l.name }"
          class="px-3 py-2 text-sm font-medium rounded-md text-gray-600 dark:text-gray-300 hover:text-gray-900 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-800"
          active-class="!text-primary-600 dark:!text-primary-400 !bg-primary-50 dark:!bg-primary-950/40"
        >
          {{ l.label }}
        </router-link>
      </nav>

      <div class="flex items-center gap-1">
        <a
          v-for="s in profile.socials.filter((s) => s.icon !== 'mail')"
          :key="s.label"
          :href="s.href"
          target="_blank"
          rel="noopener"
          :aria-label="s.label"
          class="hidden sm:inline-flex p-2 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-900 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-800"
        >
          <SocialIcon :name="s.icon" />
        </a>

        <button
          type="button"
          @click="toggle"
          class="p-2 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-900 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-800"
          :aria-label="isDark ? 'Switch to light mode' : 'Switch to dark mode'"
        >
          <svg
            v-if="isDark"
            xmlns="http://www.w3.org/2000/svg"
            viewBox="0 0 24 24"
            fill="none"
            stroke="currentColor"
            stroke-width="2"
            stroke-linecap="round"
            stroke-linejoin="round"
            class="w-5 h-5"
          >
            <circle cx="12" cy="12" r="4" />
            <path
              d="M12 2v2M12 20v2M4.93 4.93l1.41 1.41M17.66 17.66l1.41 1.41M2 12h2M20 12h2M4.93 19.07l1.41-1.41M17.66 6.34l1.41-1.41"
            />
          </svg>
          <svg
            v-else
            xmlns="http://www.w3.org/2000/svg"
            viewBox="0 0 24 24"
            fill="none"
            stroke="currentColor"
            stroke-width="2"
            stroke-linecap="round"
            stroke-linejoin="round"
            class="w-5 h-5"
          >
            <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z" />
          </svg>
        </button>

        <button
          type="button"
          class="md:hidden p-2 rounded-md text-gray-600 dark:text-gray-300 hover:bg-gray-100 dark:hover:bg-gray-800"
          @click="mobileOpen = !mobileOpen"
          aria-label="Toggle navigation"
        >
          <svg
            xmlns="http://www.w3.org/2000/svg"
            viewBox="0 0 24 24"
            fill="none"
            stroke="currentColor"
            stroke-width="2"
            stroke-linecap="round"
            stroke-linejoin="round"
            class="w-5 h-5"
          >
            <path v-if="!mobileOpen" d="M3 6h18M3 12h18M3 18h18" />
            <path v-else d="M6 6l12 12M6 18L18 6" />
          </svg>
        </button>
      </div>
    </div>

    <div v-if="mobileOpen" class="md:hidden border-t border-gray-200 dark:border-gray-800">
      <nav class="container-wide py-2 flex flex-col">
        <router-link
          v-for="l in links"
          :key="l.name"
          :to="{ name: l.name }"
          class="py-2 px-1 text-sm font-medium text-gray-600 dark:text-gray-300"
          active-class="!text-primary-600 dark:!text-primary-400"
          @click="mobileOpen = false"
        >
          {{ l.label }}
        </router-link>
      </nav>
    </div>
  </header>
</template>
