<script setup>
// Persistent branding for every slide of a deck, in the two letterbox gutters
// reveal leaves above and below the scaled slides: the deck's QR top-left and
// its project mark top-right, partner logos along the bottom with reveal's
// slide number after them. Purely decorative background dressing — not
// informative content, hence alt="".
//
// The QR points at the deck's own published URL, so it differs per
// presentation while the rest of the branding is shared. Each deck passes its
// own via the `qr` prop — by convention a `presentation.png` sitting in that
// deck's asset folder, e.g.
//   <DeckBrand :qr="asset('presentation.png')" />
// The default is the PyData deck's QR, kept only so decks predating this prop
// render exactly as before; new decks should always pass their own.
// `showLogos: false` drops the partner logos and the project mark, keeping
// only the QR chip. Decks that are not FLIP/KCL work — an invited lecture at
// another institute, say — carry no partner branding at all, but still want
// the scannable link back to the published deck. The QR keeps the top band on
// its own, so that gutter stays the same height either way and slides do not
// reflow between a branded and an unbranded deck.
const props = defineProps({
    qr: {
        type: String,
        default: `${import.meta.env.BASE_URL}static/img/logos/presentation.png`,
    },
    showLogos: { type: Boolean, default: true },
});

const base = `${import.meta.env.BASE_URL}static/img/logos/`;
const partnersLeft = ["KCL_logo.png", "aic_logo.png"];
const partnersRight = ["deepC_logo.png", "gstt_logo.png", "onelondon_logo.png"];
</script>

<template>
    <div class="deck-brand-top" aria-hidden="true">
        <span class="brand-chip brand-chip--qr">
            <img :src="props.qr" alt="" />
        </span>
        <span class="brand-mark" v-if="props.showLogos">
            <img :src="`${base}flip_logo.png`" alt="" />
        </span>
    </div>
    <div class="deck-brand-bottom" aria-hidden="true" v-if="props.showLogos">
        <span class="deck-brand-group">
            <span class="brand-chip" v-for="file in partnersLeft" :key="file">
                <img :src="`${base}${file}`" alt="" />
            </span>
        </span>
        <span class="deck-brand-group">
            <span class="brand-chip" v-for="file in partnersRight" :key="file">
                <img :src="`${base}${file}`" alt="" />
            </span>
        </span>
    </div>
</template>

<style scoped>
/* Both bands are full-width bars rather than independently-anchored corners.
   The rows used to be positioned `left: 1em` and `right: 1em` separately, so
   at deck widths where their combined intrinsic width exceeded the viewport
   they silently overlapped — the wide deepC chip landed on top of the
   left row's last chip and hid the QR code entirely. A single flex bar with
   space-between makes that collision impossible at any width.

   Rendered via RevealDeck's #chrome slot — siblings of .slides, so these
   anchor to the stable, full-size .reveal box instead of any one slide's own
   (content-height-dependent, center:true-shifted) section box. That's what
   keeps them in the same place on every slide regardless of content.

   Placement: the bottom bar used to float 2.8em up, *inside* the slide area,
   which is why tall slides ran over it and why .slide-body needed a veil.
   reveal's `margin` config (see RevealDeck.vue) subtracts a fraction of the
   viewport and centres what's left, leaving an equal letterbox gutter above
   and below the slides; both bars are sized to sit in those gutters, so slide
   content and branding can no longer occupy the same pixels. The
   viewport-relative units are what keep the two tracking each other — the
   gutters grow with the viewport, and so do these.

   The band height is set by its tallest child plus breathing room: the QR up
   top, a partner chip below. Keep them in that order if either is resized. */
.deck-brand-top,
.deck-brand-bottom {
    position: absolute;
    left: 0;
    right: 0;
    height: clamp(64px, 8.5vh, 100px);
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 0.6em;
    padding: 0 1.1em;
    pointer-events: none;
    /* Above .slides (z-index 1): these occupy reserved gutter space, so there is
     nothing underneath them to hide, and a slide can no longer drop text on
     top of a logo. */
    z-index: 2;
}
.deck-brand-top {
    margin-left: -1.2em;
    top: 0;
}
.deck-brand-bottom {
    bottom: 0;
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
    padding: 0.1em 0.4em;
    min-width: 0;
}
.brand-chip img {
    display: block;
    /* Shrinks rather than overflowing once the bar runs out of room, so a wide
     logo can never push a neighbour off-screen or under another chip. */
    height: clamp(30px, 6.2vh, 72px);
    width: auto;
    max-width: 100%;
    object-fit: contain;
}
/* The QR has the top-left corner to itself, so it can run larger than the
   partner logos — it is the one mark that has to stay scannable from the back
   of a room. It sets the band height; nothing else may exceed it. */
.brand-chip--qr img {
    height: clamp(2em, 10em, 3em);
    padding: 0.2em 0em;
    margin-top: 30px;
}
/* The deck's own mark, alone in the top-right corner. No white chip behind it:
   it is a transparent PNG that reads on both the light and dark deck grounds. */
.brand-mark {
    display: inline-flex;
    align-items: center;
    min-width: 0;
    /* Reserves the space reveal paints its slide number into (see .slide-number
     in RevealDeck.vue) so the two read as one unit: mark, then count. */
    margin-right: 3em;
}
.brand-mark img {
    display: block;
    height: clamp(32px, 6.6vh, 78px);
    width: auto;
    opacity: 0.9;
}
</style>
