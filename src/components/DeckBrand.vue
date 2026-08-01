<script setup>
// Persistent branding watermark dropped into every slide of a deck: partner
// logos anchored bottom-left/bottom-right, the deck's own project mark
// top-right. Purely decorative background dressing — not informative
// content, hence alt="".
const base = `${import.meta.env.BASE_URL}static/img/logos/`
const partnersLeft = ['KCL_logo.png', 'aic_logo.png', 'presentation.png']
const partnersRight = ['deepC_logo.png', 'gstt_logo.png', 'onelondon_logo.png']
</script>

<template>
  <div class="deck-brand-bottom deck-brand-bottom--left" aria-hidden="true">
    <span class="brand-chip" v-for="file in partnersLeft" :key="file">
      <img :src="`${base}${file}`" alt="" />
    </span>
  </div>
  <div class="deck-brand-bottom deck-brand-bottom--right" aria-hidden="true">
    <span class="brand-chip" v-for="file in partnersRight" :key="file">
      <img :src="`${base}${file}`" alt="" />
    </span>
  </div>
  <div class="deck-brand-top-right" aria-hidden="true">
    <img :src="`${base}flip_logo.png`" alt="" />
  </div>
</template>

<style scoped>
.deck-brand-bottom {
  /* Rendered via RevealDeck's #chrome slot — a sibling of .slides, so this
     anchors to the stable, full-size .reveal box instead of any one slide's
     own (content-height-dependent, center:true-shifted) section box. That's
     what keeps it in the same place on every slide regardless of content. */
  position: absolute;
  bottom: 2.8em;
  display: flex;
  align-items: center;
  gap: 0.5em;
  pointer-events: none;
}
.deck-brand-bottom--left {
  left: 1em;
}
.deck-brand-bottom--right {
  right: 1em;
}
.brand-chip {
  display: inline-flex;
  align-items: center;
  background: #ffffff;
  border-radius: 5px;
  padding: 0.2em 0.55em;
}
.brand-chip img {
  display: block;
  height: 2.6em;
  width: auto;
}
.deck-brand-top-right {
  position: absolute;
  top: 0.6em;
  right: 0.8em;
  pointer-events: none;
}
.deck-brand-top-right img {
  display: block;
  height: 3em;
  width: auto;
  opacity: 0.9;
  filter: drop-shadow(0 2px 6px rgba(0, 0, 0, 0.3));
}
</style>
