<script setup>
// Persistent branding watermark dropped into every slide of a deck: partner
// logos anchored bottom-left/bottom-right, the deck's own project mark
// top-right. Purely decorative background dressing — not informative
// content, hence alt="".
//
// The last chip on the left is a QR code pointing at the deck's own published
// URL, so it differs per presentation while the rest of the branding is shared.
// Each deck passes its own via the `qr` prop — by convention a `presentation.png`
// sitting in that deck's asset folder, e.g.
//   <DeckBrand :qr="asset('presentation.png')" />
// The default is the PyData deck's QR, kept only so decks predating this prop
// render exactly as before; new decks should always pass their own.
//
// `showLogos` hides the partner/project logos (KCL, AI Centre, GSTT, deepC,
// One London, and the FLIP mark) while keeping the deck's own QR. For decks
// that are not part of the KCL/FLIP work — e.g. the IAA-SO astronomy lecture.
const props = defineProps({
  qr: { type: String, default: `${import.meta.env.BASE_URL}static/img/logos/presentation.png` },
  showLogos: { type: Boolean, default: true }
})

const base = `${import.meta.env.BASE_URL}static/img/logos/`
const partnersLeft = ['KCL_logo.png', 'aic_logo.png']
const partnersRight = ['deepC_logo.png', 'gstt_logo.png', 'onelondon_logo.png']
</script>

<template>
  <div class="deck-brand-bottom" aria-hidden="true">
    <span class="deck-brand-group">
      <template v-if="showLogos">
        <span class="brand-chip" v-for="file in partnersLeft" :key="file">
          <img :src="`${base}${file}`" alt="" />
        </span>
      </template>
      <span class="brand-chip brand-chip--qr">
        <img :src="props.qr" alt="" />
      </span>
    </span>
    <span v-if="showLogos" class="deck-brand-group">
      <span class="brand-chip" v-for="file in partnersRight" :key="file">
        <img :src="`${base}${file}`" alt="" />
      </span>
    </span>
  </div>
  <div v-if="showLogos" class="deck-brand-top-right" aria-hidden="true">
    <img :src="`${base}flip_logo.png`" alt="" />
  </div>
</template>

<style scoped>
/* One full-width bar rather than two independently-anchored rows.
   The rows used to be positioned `left: 1em` and `right: 1em` separately, so
   at deck widths where their combined intrinsic width exceeded the viewport
   they silently overlapped — the wide deepC chip landed on top of the
   left row's last chip and hid the QR code entirely. A single flex bar with
   space-between makes that collision impossible at any width. */
.deck-brand-bottom {
  /* Rendered via RevealDeck's #chrome slot — a sibling of .slides, so this
     anchors to the stable, full-size .reveal box instead of any one slide's
     own (content-height-dependent, center:true-shifted) section box. That's
     what keeps it in the same place on every slide regardless of content. */
  position: absolute;
  bottom: 2.8em;
  left: 1em;
  right: 1em;
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 0.6em;
  pointer-events: none;
  /* Sits between reveal's two layers: above .backgrounds (z-index auto), so an
     opaque background such as a background video can't hide the branding, but
     below .slides (also z-index 1, and later in DOM order, so it wins the
     tie-break) — slide text still paints over the logos where the two overlap,
     and the .slide-body veil is what keeps that text readable. */
  z-index: 1;
}
.deck-brand-group {
  display: flex;
  align-items: center;
  gap: 0.5em;
  min-width: 0;
}
.brand-chip {
  display: inline-flex;
  align-items: center;
  background: #ffffff;
  border-radius: 5px;
  padding: 0.2em 0.55em;
  min-width: 0;
}
.brand-chip img {
  display: block;
  /* Shrinks rather than overflowing once the bar runs out of room, so a wide
     logo can never push a neighbour off-screen or under another chip. */
  height: 2.3em;
  width: auto;
  max-width: 100%;
  object-fit: contain;
}
/* The QR is the one chip that has to stay scannable from the back of a room,
   so it keeps its size while the partner logos absorb any shortfall. */
.brand-chip--qr img {
  height: 2.6em;
  flex: none;
}
.deck-brand-top-right {
  position: absolute;
  top: 0.6em;
  right: 0.8em;
  pointer-events: none;
  z-index: 1;
}
.deck-brand-top-right img {
  display: block;
  height: 3em;
  width: auto;
  opacity: 0.9;
  filter: drop-shadow(0 2px 6px rgba(0, 0, 0, 0.3));
}
</style>
