<script setup>
// Present mode: one deck, no router, no site chrome. Reached at /?deck=<id>.
//
// Why it exists: reveal's speaker view (press S) addresses slides through the
// URL hash — it rebuilds the deck URL as "<path>?receiver&…#/12/0" for its two
// preview iframes. Everywhere else in this SPA the hash belongs to vue-router,
// so those iframes would load the home page instead of the deck, never post
// their "ready" event, and the speaker view would sit on "Loading speaker
// view…" for ever with no notes and no sync. Keying the deck off the query
// string instead leaves the hash free for reveal, which is all it takes.
//
// The normal /#/presentations/:id route is unaffected and still the one to
// share; this is the one to present from.
import { computed, defineAsyncComponent, onMounted, shallowRef, watchEffect } from 'vue'
import { getPresentation } from '@/content/presentations'
import { useDarkMode } from '@/composables/useDarkMode'

const props = defineProps({
  id: { type: String, required: true }
})

// The site header — and with it the theme toggle — never mounts here, so apply
// the stored light/dark preference directly.
useDarkMode()

const meta = computed(() => getPresentation(props.id))
const deckComponent = shallowRef(null)

watchEffect(() => {
  const p = meta.value
  if (p) deckComponent.value = defineAsyncComponent(p.deck)
})

onMounted(() => {
  document.title = meta.value ? `${meta.value.title} — presenting` : 'Deck not found'
})
</script>

<template>
  <!-- reveal measures its own container, so the deck needs a viewport-sized box
       here exactly as it gets one on the router route. Without it .reveal is
       0px tall, reveal floors the scale at minScale, and the slide paints as a
       sliver — which also showed up as blank previews in the speaker view. -->
  <div v-if="deckComponent" class="fixed inset-0 bg-white dark:bg-dracula-bg">
    <component :is="deckComponent" />
  </div>
  <div v-else class="flex min-h-screen items-center justify-center p-8 text-center">
    <div>
      <h1 class="text-2xl font-bold text-gray-900 dark:text-white">No deck with id "{{ id }}"</h1>
      <p class="mt-2 text-gray-600 dark:text-gray-400">
        Present mode takes a presentation id, e.g.
        <code>/?deck=flip-reusable-fl-platform-2026</code>.
      </p>
      <a class="link mt-4 inline-block" href="./#/presentations">Back to the presentations</a>
    </div>
  </div>
</template>
