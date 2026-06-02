<script setup>
import { computed, defineAsyncComponent, ref, shallowRef, watchEffect, onMounted, onBeforeUnmount } from 'vue'
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

// --- Presenter chrome: auto-hide on mouse stillness while fullscreen ---------
// In fullscreen, the Exit button and reveal's nav arrows fade out after a short
// idle period and reappear the moment the mouse moves. Outside fullscreen the
// controls always stay visible so there's always an obvious way out.
const isFullscreen = ref(false)
const idle = ref(false)
let idleTimer = null

const hideChrome = computed(() => isFullscreen.value && idle.value)

// --- Laser pointer (toggle with "L") -----------------------------------------
// The dot only shows while the mouse is moving; it vanishes ~0.3s after the
// mouse goes still, like a real laser you stop pointing.
const laserOn = ref(false)
const laserMoving = ref(false)
const laser = ref({ x: 0, y: 0 })
let laserTimer = null

const showLaser = computed(() => laserOn.value && laserMoving.value)

function wake(e) {
  if (e && e.type === 'mousemove') {
    laser.value = { x: e.clientX, y: e.clientY }
    laserMoving.value = true
    clearTimeout(laserTimer)
    laserTimer = setTimeout(() => {
      laserMoving.value = false
    }, 300)
  }
  idle.value = false
  clearTimeout(idleTimer)
  idleTimer = setTimeout(() => {
    idle.value = true
  }, 500)
}

function onFullscreenChange() {
  isFullscreen.value = !!document.fullscreenElement
  wake()
}

function onKeydown(e) {
  // Don't steal typing in inputs; the deck has none, but be safe.
  if (e.key === 'l' || e.key === 'L') laserOn.value = !laserOn.value
}

onMounted(() => {
  window.addEventListener('mousemove', wake)
  window.addEventListener('mousedown', wake)
  window.addEventListener('keydown', onKeydown)
  document.addEventListener('fullscreenchange', onFullscreenChange)
  isFullscreen.value = !!document.fullscreenElement
  wake()
})

onBeforeUnmount(() => {
  window.removeEventListener('mousemove', wake)
  window.removeEventListener('mousedown', wake)
  window.removeEventListener('keydown', onKeydown)
  document.removeEventListener('fullscreenchange', onFullscreenChange)
  clearTimeout(idleTimer)
  clearTimeout(laserTimer)
})
</script>

<template>
  <!-- Fullscreen overlay so reveal.js controls the whole viewport. Background
       follows the site theme (white in light mode, Dracula bg in dark). -->
  <div
    class="fixed inset-0 z-50 bg-white dark:bg-dracula-bg deck-overlay"
    :class="{ 'chrome-hidden': hideChrome, 'laser-on': laserOn }"
  >
    <component :is="deckComponent" v-if="deckComponent" />

    <!-- Presenter controls (bottom-left). Auto-hide with .chrome-hidden. -->
    <div class="deck-chrome absolute bottom-4 left-4 z-[60] flex items-center gap-2">
      <router-link
        :to="{ name: 'presentations' }"
        class="inline-flex items-center gap-2 rounded-full bg-black/40 hover:bg-black/60 backdrop-blur px-3.5 py-2 text-xs font-medium text-white ring-1 ring-white/20 transition"
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

      <button
        type="button"
        @click="laserOn = !laserOn"
        :class="laserOn ? 'bg-red-500/80 ring-red-300/40' : 'bg-black/40 hover:bg-black/60 ring-white/20'"
        class="inline-flex items-center gap-2 rounded-full backdrop-blur px-3.5 py-2 text-xs font-medium text-white ring-1 transition"
        :title="laserOn ? 'Laser pointer on (L)' : 'Laser pointer (L)'"
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
          <circle cx="12" cy="12" r="3" />
          <path d="M12 2v3M12 19v3M2 12h3M19 12h3" />
        </svg>
        Pointer
      </button>
    </div>

    <!-- Laser dot — follows the cursor while pointer mode is on, and vanishes
         ~0.1s after the mouse goes still. -->
    <div
      v-if="showLaser"
      class="laser-dot"
      :style="{ left: laser.x + 'px', top: laser.y + 'px' }"
    ></div>
  </div>
</template>

<style scoped>
/* Smoothly fade the presenter chrome and reveal's nav arrows. */
.deck-chrome,
.deck-overlay :deep(.reveal .controls) {
  transition: opacity 0.3s ease;
}

/* Idle while fullscreen: hide cursor, Exit/Pointer buttons and the nav arrows. */
.deck-overlay.chrome-hidden {
  cursor: none;
}
.deck-overlay.chrome-hidden .deck-chrome {
  opacity: 0;
  pointer-events: none;
}
.deck-overlay.chrome-hidden :deep(.reveal .controls) {
  opacity: 0;
  pointer-events: none;
}

/* Laser pointer mode: hide the system cursor; the dot shows where it is. */
.deck-overlay.laser-on {
  cursor: none;
}
.laser-dot {
  position: fixed;
  width: 20px;
  height: 20px;
  margin: -10px 0 0 -10px;
  border-radius: 50%;
  background: radial-gradient(
    circle,
    rgba(255, 90, 90, 0.95) 0%,
    rgba(255, 0, 0, 0.55) 45%,
    rgba(255, 0, 0, 0) 72%
  );
  box-shadow: 0 0 14px 5px rgba(255, 0, 0, 0.45);
  pointer-events: none;
  z-index: 70;
}
</style>
