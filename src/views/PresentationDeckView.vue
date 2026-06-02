<script setup>
import { computed, defineAsyncComponent, shallowRef, watchEffect } from 'vue'
import { useRoute, useRouter } from 'vue-router'
import { getPresentation } from '@/content/presentations'

const route = useRoute()
const router = useRouter()

const meta = computed(() => getPresentation(route.params.id))
const deckComponent = shallowRef(null)

// Resolve the lazy deck import for the requested id. Redirect to the gallery if
// the id is unknown.
watchEffect(() => {
  const p = meta.value
  if (!p) {
    router.replace({ name: 'presentations' })
    return
  }
  deckComponent.value = defineAsyncComponent(p.deck)
})
</script>

<template>
  <!-- Fullscreen overlay so reveal.js controls the whole viewport. Background
       follows the site theme (white in light mode, Dracula bg in dark). -->
  <div class="fixed inset-0 z-50 bg-white dark:bg-dracula-bg">
    <component :is="deckComponent" v-if="deckComponent" />

    <!-- Exit back to the gallery -->
    <router-link
      :to="{ name: 'presentations' }"
      class="absolute top-4 left-4 z-[60] inline-flex items-center gap-2 rounded-full bg-black/40 hover:bg-black/60 backdrop-blur px-3.5 py-2 text-xs font-medium text-white ring-1 ring-white/20 transition"
      title="Back to presentations"
    >
      <svg
        xmlns="http://www.w3.org/2000/svg"
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        stroke-width="2"
        stroke-linecap="round"
        stroke-linejoin="round"
        class="w-4 h-4"
      >
        <path d="M19 12H5M12 19l-7-7 7-7" />
      </svg>
      Exit
    </router-link>
  </div>
</template>
