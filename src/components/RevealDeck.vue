<script setup>
import Reveal from 'reveal.js'
import 'reveal.js/reveal.css'
import { nextTick, onBeforeUnmount, onMounted, ref } from 'vue'

// A thin Vue wrapper around reveal.js. Slot content must be a sequence of
// <section> elements (the slides). The deck is initialised on mount and torn
// down on unmount so navigating away leaves no global key/resize listeners.
const props = defineProps({
  options: { type: Object, default: () => ({}) }
})

const root = ref(null)
let deck = null

// Present mode (/?deck=<id>) mounts a deck with no router in the page, so reveal
// can own the hash there — and it has to, because the speaker view addresses
// slides by hash. On the /#/presentations/:id route the hash is vue-router's.
const presentMode = new URLSearchParams(window.location.search).has('deck')

onMounted(async () => {
  await nextTick()

  // Speaker view: press S for a second window carrying the notes, a clock and
  // the next slide, kept in sync with this one over postMessage — so the deck
  // can be full-screen on the projector while the notes stay on the laptop, and
  // either window can advance the talk. reveal 6 inlines the speaker-view markup
  // and document.write()s it into a popup, so it needs no plugin file on disk and
  // works as-is through Vite; it does need popups allowed for the origin.
  //
  // Loaded only in present mode, for two reasons: on the router-driven route the
  // speaker view cannot reach the deck (its preview iframes address slides by
  // hash, which is vue-router's there), so S would open a window that never
  // connects — and the plugin bundles a markdown parser, which is dead weight on
  // the route visitors actually land on.
  const plugins = presentMode ? [(await import('reveal.js/plugin/notes')).default] : []

  deck = new Reveal(root.value, {
    plugins,
    // On the router-driven route reveal's own hash navigation would fight with
    // vue-router, so it stays off there and the deck is driven by the keyboard.
    hash: presentMode,
    embedded: false,
    controls: true,
    progress: true,
    slideNumber: 'c/t',
    overview: true,
    center: false,
    transition: 'slide',
    transitionSpeed: 'default',
    width: 1280,
    height: 720,
    margin: 0.055,
    minScale: 0.2,
    maxScale: 1.8,
    ...props.options
  })
  await deck.initialize()
})

onBeforeUnmount(() => {
  try {
    deck?.destroy()
  } catch (e) {
    // reveal occasionally throws on destroy if init was interrupted — ignore.
  }
  deck = null
})
</script>

<template>
  <div ref="root" class="reveal deck-theme">
    <!-- Optional persistent overlay (e.g. branding watermarks): a sibling of
         .slides rather than nested in a section, so it anchors to the stable,
         full-size .reveal box instead of any one slide's own (content-height-
         dependent, center:true-shifted) box. Renders nothing unless a deck
         fills the slot.

         Deliberately placed BEFORE .slides. reveal gives .slides z-index 1 and
         leaves .backgrounds on auto, so chrome at z-index 1 paints above the
         backgrounds (beating an opaque background video) while losing the
         same-z tie-break to .slides, which comes later in DOM order — meaning
         slide text still paints over the logos. Any higher z-index would put
         the chrome above the slide content and hide it. -->
    <slot name="chrome" />
    <div class="slides">
      <slot />
    </div>
  </div>
</template>

<!--
  Custom reveal.js theme using the site's Dracula palette + fonts, namespaced
  under .deck-theme so it never leaks into the rest of the SPA.
-->
<style>
.reveal.deck-theme {
  /* Mode-independent typography */
  --r-main-font: 'Inter', ui-sans-serif, system-ui, sans-serif;
  --r-heading-font: 'Bai Jamjuree', 'Inter', sans-serif;
  --r-code-font: 'JetBrainsMono', ui-monospace, monospace;
  font-family: var(--r-main-font);
  /* Base type scale: reveal core defaults to 20pt; bump +20% for legibility
     when projected. All slide styles are em-relative, so this scales uniformly. */
  font-size: 24pt;
  color: var(--r-main-color);
}

/* DARK theme — follows the site's .dark class on <html> (Dracula palette) */
html.dark .reveal.deck-theme {
  --r-background-color: #282a36;
  --r-main-color: #f8f8f2;
  --r-heading-color: #f8f8f2;
  --r-link-color: #8be9fd;
  --r-link-color-hover: #bd93f9;
  --accent: #bd93f9;
  --accent-cyan: #8be9fd;
  --accent-pink: #ff79c6;
  --accent-green: #50fa7b;
  --accent-orange: #ffb86c;
  --surface: #21222c;
  --line: #44475a;
  --comment: #6272a4;
  --grad-a: #bd93f9;
  --grad-b: #8be9fd;
  --title-sub: #f8f8f2;
  --fig-bg: #f8f8f2;
  --body-veil: rgba(40, 42, 54, 0.82);
  --fig-border: transparent;
  --fig-shadow: 0 10px 30px rgba(0, 0, 0, 0.35);
  --flip-bg: rgba(139, 233, 253, 0.08);
  --flip-border: #8be9fd;
  --flip-tag-bg: #8be9fd;
  --flip-tag-fg: #282a36;
  --pill-bg: rgba(189, 147, 249, 0.14);
  --pill-color: #f8f8f2;
  --req-bg: rgba(189, 147, 249, 0.16);
  --fla3-col-bg: rgba(189, 147, 249, 0.08);
  --flip-col-bg: rgba(139, 233, 253, 0.07);
  --code-color: #50fa7b;
  background:
    linear-gradient(to bottom, #282A3600 0%, #282a36 100%),
}

/* LIGHT theme — default for presenting; better for holding attention */
html:not(.dark) .reveal.deck-theme {
  --r-background-color: #ffffff;
  --r-main-color: #21222c;
  --r-heading-color: #1a1b26;
  --r-link-color: #0e7490;
  --r-link-color-hover: #7c3aed;
  --accent: #7c3aed;
  --accent-cyan: #0e7490;
  --accent-pink: #db2777;
  --accent-green: #15803d;
  --accent-orange: #c2410c;
  --surface: #f4f3fb;
  --line: #d9d8e6;
  --comment: #64748b;
  --grad-a: #7c3aed;
  --grad-b: #0891b2;
  --title-sub: #21222c;
  --fig-bg: #ffffff;
  --body-veil: rgba(255, 255, 255, 0.85);
  --fig-border: #e2e2ec;
  --fig-shadow: 0 8px 24px rgba(24, 18, 48, 0.12);
  --flip-bg: rgba(8, 145, 178, 0.07);
  --flip-border: #0891b2;
  --flip-tag-bg: #0e7490;
  --flip-tag-fg: #ffffff;
  --pill-bg: rgba(124, 58, 237, 0.1);
  --pill-color: #3b2a63;
  --req-bg: rgba(124, 58, 237, 0.1);
  --fla3-col-bg: rgba(124, 58, 237, 0.07);
  --flip-col-bg: rgba(8, 145, 178, 0.07);
  --code-color: #b91c1c;
  background:
    linear-gradient(to bottom, #DF428B35 0%, #DF428B20 10%, #ffffff00 60%);
}

.reveal.deck-theme .slides {
  text-align: left;
}

.reveal.deck-theme h1,
.reveal.deck-theme h2,
.reveal.deck-theme h3,
.reveal.deck-theme h4 {
  font-family: var(--r-heading-font);
  font-weight: 700;
  letter-spacing: -0.01em;
  text-transform: none;
  line-height: 1.1;
  margin-bottom: 0.4em;
}

.reveal.deck-theme h1 { font-size: 2.1em; }
.reveal.deck-theme h2 { font-size: 1.5em; }
.reveal.deck-theme h3 { font-size: 1.15em; color: var(--accent); }

.reveal.deck-theme .eyebrow {
  font-family: var(--r-code-font);
  text-transform: uppercase;
  letter-spacing: 0.18em;
  font-size: 0.42em;
  color: var(--accent);
  margin-bottom: 0.6em;
}

.reveal.deck-theme p,
.reveal.deck-theme li {
  line-height: 1.4;
}

.reveal.deck-theme strong { color: var(--accent); font-weight: 700; }
.reveal.deck-theme em { color: var(--accent-cyan); font-style: normal; }

.reveal.deck-theme .accent { color: var(--accent); }
.reveal.deck-theme .cyan { color: var(--accent-cyan); }
.reveal.deck-theme .pink { color: var(--accent-pink); }
.reveal.deck-theme .green { color: var(--accent-green); }
.reveal.deck-theme .orange { color: var(--accent-orange); }
.reveal.deck-theme .muted { color: var(--comment); }
.reveal.deck-theme .small { font-size: 0.7em; }
.reveal.deck-theme .medium { font-size: 0.85em; }
.reveal.deck-theme .center { text-align: center; }

/* Title slide */
.reveal.deck-theme .title-slide h1 {
  font-size: 2.0em;
  background: linear-gradient(120deg, var(--grad-a), var(--grad-b));
  -webkit-background-clip: text;
  background-clip: text;
  -webkit-text-fill-color: transparent;
}
.reveal.deck-theme .title-slide .subtitle {
  font-size: 1.05em;
  color: var(--title-sub);
  margin-top: 0.2em;
}
.reveal.deck-theme .title-slide .byline {
  margin-top: 0.2em;
  color: var(--comment);
  font-size: 0.7em;
}
.reveal.deck-theme .title-slide .venue-note {
  display: inline-block;
  margin-top: 0.2em;
  font-family: var(--r-code-font);
  font-size: 0.6em;
  letter-spacing: 0.04em;
  color: var(--accent);
  background: var(--req-bg);
  border: 1px solid var(--line);
  border-radius: 999px;
  padding: 0.3em 1em;
}

/* Cards / columns */
.reveal.deck-theme .cols {
  display: grid;
  grid-template-columns: repeat(var(--n, 2), 1fr);
  gap: 0.8em;
  align-items: stretch;
}
.reveal.deck-theme .panel {
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: 14px;
  padding: 0.7em 0.9em;
}
.reveal.deck-theme .panel h3 { margin-top: 0; }
/* Compact row: panels that are a supporting inventory rather than the slide's
   argument — a label, one short line, its pills — so the figure above them can
   have the vertical space instead. */
.reveal.deck-theme .cols.compact .panel { padding: 0.45em 0.7em; }
.reveal.deck-theme .cols.compact .panel h3 { font-size: 0.85em; margin-bottom: 0.1em; }
.reveal.deck-theme .cols.compact .panel p { margin: 0 0 0.2em; }
.reveal.deck-theme .cols.compact .panel .pill { margin-top: 0; margin-bottom: 0; }
.reveal.deck-theme .panel.fla3 { border-top: 3px solid var(--accent); }
.reveal.deck-theme .panel.flip { border-top: 3px solid var(--accent-cyan); }

/* Pills */
.reveal.deck-theme .pill {
  display: inline-block;
  background: var(--pill-bg);
  border: 1px solid var(--line);
  color: var(--pill-color);
  border-radius: 999px;
  padding: 0.12em 0.7em;
  font-size: 0.6em;
  font-family: var(--r-code-font);
  margin: 0.15em 0.2em 0.15em 0;
}
/* The pill that is us, in a list of other people's work: same shape and size,
   so it reads as one of the group rather than a ranking, but findable at a
   glance from the back of the room. Accent border and text rather than an
   accent fill, which keeps it legible in both themes. */
.reveal.deck-theme .pill.ours {
  border-color: var(--accent);
  color: var(--accent);
  font-weight: 700;
}
/* A dot marking pills that carry a second, orthogonal fact — on a roll-call of
   other projects, "someone from this one is in the room". Kept as a ::before
   marker rather than another colour so it composes with .ours instead of
   fighting it over the same two properties. .here-dot is the same mark, for
   naming the convention in a legend line. */
.reveal.deck-theme .pill.here::before,
.reveal.deck-theme .here-dot {
  content: '';
  display: inline-block;
  width: 0.5em;
  height: 0.5em;
  border-radius: 50%;
  background: var(--accent-cyan);
  margin-right: 0.4em;
}
.reveal.deck-theme .here-dot { width: 0.42em; height: 0.42em; }

/* Tables */
.reveal.deck-theme table {
  font-size: 0.62em;
  border-collapse: collapse;
  margin: 0.3em 0;
  width: 100%;
}
.reveal.deck-theme table th,
.reveal.deck-theme table td {
  border: 1px solid var(--line);
  padding: 0.4em 0.6em;
  text-align: left;
  vertical-align: top;
}
.reveal.deck-theme table th {
  background: var(--surface);
  font-family: var(--r-heading-font);
}
.reveal.deck-theme table td.fla3-col { background: var(--fla3-col-bg); }
.reveal.deck-theme table td.flip-col { background: var(--flip-col-bg); }

/* Code */
.reveal.deck-theme code {
  font-family: var(--r-code-font);
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: 6px;
  padding: 0.05em 0.35em;
  font-size: 0.85em;
  color: var(--code-color);
}
.reveal.deck-theme pre code {
  display: block;
  padding: 0.8em 1em;
  line-height: 1.4;
}

/* Big stat */
.reveal.deck-theme .stat {
  font-family: var(--r-heading-font);
  font-size: 2.2em;
  font-weight: 700;
  color: var(--accent);
  line-height: 1;
}
.reveal.deck-theme .stat-label { font-size: 0.55em; color: var(--comment); }

/* Controls / progress tint */
.reveal.deck-theme .controls { color: var(--accent); }
.reveal.deck-theme .progress { color: var(--accent); height: 4px; }
.reveal.deck-theme .slide-number { background: transparent; color: var(--comment); }

/* Lists */
.reveal.deck-theme ul { margin-left: 1em; }
.reveal.deck-theme li { margin-bottom: 0.35em; }
.reveal.deck-theme .checklist { list-style: none; margin-left: 0; }
.reveal.deck-theme .checklist li::before {
  content: '\2713';
  color: var(--accent-green);
  font-weight: 700;
  margin-right: 0.5em;
}
.reveal.deck-theme .crosslist { list-style: none; margin-left: 0; }
.reveal.deck-theme .crosslist li::before {
  content: '\2717';
  color: var(--accent-pink);
  font-weight: 700;
  margin-right: 0.5em;
}
.reveal.deck-theme .dotlist { list-style: none; margin-left: 0; }
.reveal.deck-theme .dotlist li::before {
  content: '\2022';
  color: var(--accent);
  font-weight: 700;
  margin-right: 0.5em;
}

/* Compact reference glossary (appendix) — dense multi-column term list, denser
   than content slides on purpose since it is scanned, not presented line-by-line. */
.reveal.deck-theme .glossary {
  column-count: 3;
  column-gap: 1.5em;
  font-size: 0.5em;
  line-height: 1.3;
  text-align: left;
}
.reveal.deck-theme .glossary .gl-group {
  break-inside: avoid;
  margin-bottom: 0.7em;
}
.reveal.deck-theme .glossary h4 {
  color: var(--accent);
  font-family: var(--r-code-font);
  font-size: 1em;
  text-transform: uppercase;
  letter-spacing: 0.07em;
  margin: 0 0 0.3em;
}
.reveal.deck-theme .glossary p { margin: 0 0 0.18em; }
.reveal.deck-theme .glossary b { color: var(--accent-cyan); font-weight: 700; }

/* Horizontal "worry scale" — gradient track with 1–5 ticks under three meme cards. */
.reveal.deck-theme .scale-bar {
  height: 0.4em;
  border-radius: 999px;
  margin: 0.2em 0 0.15em;
  background: linear-gradient(90deg, var(--accent-green) 0%, var(--accent-orange) 50%, var(--accent-pink) 100%);
}
.reveal.deck-theme .scale-ticks {
  display: flex;
  justify-content: space-between;
  font-family: var(--r-code-font);
  font-size: 0.55em;
  color: var(--comment);
}
.reveal.deck-theme .worry-num {
  font-family: var(--r-heading-font);
  font-weight: 700;
  font-size: 1.5em;
  line-height: 1;
}

/* Everything on a slide that isn't the eyebrow or the heading lives in one of
   these. The persistent DeckBrand logo bar is anchored to the .reveal box, so
   tall slides can run over it — white logo chips behind body text made both
   unreadable. A semi-transparent veil in the deck's own background colour keeps
   the text legible wherever it lands, while still letting the branding show
   through. Backdrop blur is a progressive enhancement: browsers without it
   simply get the flat translucent panel, which is already enough. */
.reveal.deck-theme .slide-body {
  background: var(--body-veil);
  backdrop-filter: blur(4px);
  -webkit-backdrop-filter: blur(4px);
  border-radius: 14px;
  padding: 0.5em 0.7em;
  position: relative;
  z-index: 1;
}
/* Centred slides (title/divider) read better without a hard-edged panel, so the
   veil there is a soft pad rather than a card. */
.reveal.deck-theme .title-slide .slide-body {
  background: none;
  backdrop-filter: none;
  -webkit-backdrop-filter: none;
  padding: 0;
}
/* fit-content + auto margins, NOT inline-block: the veil hugs the text but each
   child keeps its own line, so a byline and a pill don't end up side by side. */
.reveal.deck-theme .title-slide .slide-body > * {
  width: fit-content;
  max-width: 100%;
  margin-left: auto;
  margin-right: auto;
  background: var(--body-veil);
  backdrop-filter: blur(4px);
  -webkit-backdrop-filter: blur(4px);
  border-radius: 12px;
  padding: 0.15em 0.6em;
}
/* .venue-note is already a pill with its own background and padding — leave it
   alone rather than stacking a second veil behind it. */
.reveal.deck-theme .title-slide .slide-body > .venue-note {
  background: var(--req-bg);
  padding: 0.3em 1em;
}

/* A slide whose backdrop is a full-bleed video (reveal's own
   data-background-video layer). The FLIP globe artwork is light, so this slide
   flips to ink-on-white in BOTH themes rather than inheriting the dark palette —
   otherwise pale text would sit on a near-white globe. Content rides on a soft
   card so the dotted globe never eats a letter. */
/* Poster behind the background video. reveal offers no poster option for
   data-background-video, and setting data-background-image instead would win the
   if/else in its background loader and suppress the video entirely — so the
   still is painted on the background layer underneath, showing through until the
   video paints, or permanently if the browser blocks autoplay. */
.reveal.deck-theme .slide-background.video-hero .slide-background-content {
  background-image: var(--hero-poster);
  background-size: cover;
  background-position: center;
}
.reveal.deck-theme section.video-hero {
  --r-main-color: #21222c;
  --r-heading-color: #1a1b26;
  --accent: #a21caf;
  --accent-cyan: #0e7490;
  --comment: #52525b;
  --grad-a: #7c3aed;
  --grad-b: #db2777;
  --req-bg: rgba(124, 58, 237, 0.12);
  --line: #cfcde0;
  --title-sub: #21222c;
  --body-veil: transparent;
  color: var(--r-main-color);
  background: rgba(255, 255, 255, 0.62);
  backdrop-filter: blur(3px);
  -webkit-backdrop-filter: blur(3px);
  border-radius: 18px;
  padding: 0.6em 1em 0.8em;
}
/* The section is already the card — don't stack a second veil inside it. */
.reveal.deck-theme section.video-hero .slide-body > * {
  background: none;
  backdrop-filter: none;
  -webkit-backdrop-filter: none;
}
.reveal.deck-theme section.video-hero .slide-body > .venue-note {
  background: var(--req-bg);
}

/* Figures from the paper — framed on a light card. In dark mode the card lifts
   the matplotlib white backgrounds off the dark deck; in light mode a subtle
   border separates the (white) figure from the (white) slide. */
.reveal.deck-theme .figure {
  background: var(--fig-bg);
  border: 1px solid var(--fig-border);
  border-radius: 12px;
  padding: 0.5em;
  display: inline-block;
  box-shadow: var(--fig-shadow);
}
.reveal.deck-theme .figure img,
.reveal.deck-theme .figure video {
  display: block;
  margin: 0;
  border-radius: 6px;
  max-height: 100%;
}
.reveal.deck-theme figcaption,
.reveal.deck-theme .caption {
  font-size: 0.4em;
  color: var(--comment);
  margin-top: 0.5em;
  font-family: var(--r-code-font);
}

/* Woven-in FLIP contrast aside — appears on slides where FLA3's choice has a
   direct FLIP counterpart. Always cyan-keyed so the audience learns the cue. */
.reveal.deck-theme .flip-note {
  border: 1px solid var(--flip-border);
  border-left: 5px solid var(--flip-border);
  background: var(--flip-bg);
  border-radius: 8px;
  padding: 0.75em 0.95em;
  font-size: 0.8em;
  line-height: 1.45;
  margin-top: 0.95em;
}
.reveal.deck-theme .flip-note .flip-tag {
  display: inline-block;
  font-family: var(--r-code-font);
  font-weight: 700;
  letter-spacing: 0.08em;
  color: var(--flip-tag-fg);
  background: var(--flip-tag-bg);
  border-radius: 5px;
  padding: 0.1em 0.55em;
  margin-right: 0.55em;
  font-size: 0.92em;
}

/* Compact FLIP comparison line living *inside* a panel or list item — shows,
   per item, whether FLIP covers the same ground. Cyan-keyed like the
   block-level note, separated from the body by a dashed rule. */
.reveal.deck-theme .flip-inline {
  margin-top: 0.7em;
  padding-top: 0.55em;
  border-top: 1px dashed var(--flip-border);
  font-size: 0.62em;
  line-height: 1.4;
}
/* Tighter inside the contributions list — 5 rows, keep it on-slide. */
.reveal.deck-theme ol.contribs .flip-inline {
  margin-top: 0.35em;
  padding-top: 0.3em;
  font-size: 0.8em;
}
/* Requirements list (R1–R5): five items each gaining a FLIP verdict, so pull
   the list type down slightly and tighten every gap to keep it on-slide. */
.reveal.deck-theme ul.reqs { font-size: 0.86em; }
.reveal.deck-theme ul.reqs > li { margin-bottom: 0.2em; line-height: 1.3; }
.reveal.deck-theme ul.reqs .flip-inline {
  margin-top: 0.22em;
  padding-top: 0.22em;
  font-size: 0.68em;
}
.reveal.deck-theme .flip-inline .flip-tag {
  display: inline-block;
  font-family: var(--r-code-font);
  font-weight: 700;
  letter-spacing: 0.08em;
  color: var(--flip-tag-fg);
  background: var(--flip-tag-bg);
  border-radius: 4px;
  padding: 0.04em 0.4em;
  margin-right: 0.45em;
  font-size: 0.95em;
}
.reveal.deck-theme .cov {
  font-weight: 700;
  margin-right: 0.3em;
}
.reveal.deck-theme .cov.yes { color: var(--accent-green); }
.reveal.deck-theme .cov.partial { color: var(--accent-orange); }
.reveal.deck-theme .cov.no { color: var(--accent-pink); }

/* Numbered contribution list with accent number badges */
.reveal.deck-theme ol.contribs {
  list-style: none;
  margin-left: 0;
  counter-reset: c;
}
.reveal.deck-theme ol.contribs > li {
  counter-increment: c;
  position: relative;
  padding-left: 2em;
  margin-bottom: 0.4em;
}
/* Five-item version of the list: the default gap tips it over the slide box, and
   losing a third of the gap is cheaper than losing type size. */
.reveal.deck-theme ol.contribs.tight > li { margin-bottom: 0.22em; }
.reveal.deck-theme ol.contribs > li::before {
  content: counter(c);
  position: absolute;
  left: 0;
  top: 0.05em;
  font-family: var(--r-heading-font);
  font-weight: 700;
  color: var(--accent);
  background: var(--req-bg);
  border: 1px solid var(--line);
  border-radius: 7px;
  width: 1.4em;
  height: 1.4em;
  display: inline-flex;
  align-items: center;
  justify-content: center;
  font-size: 0.9em;
}

/* Requirement chips R1–R5 */
.reveal.deck-theme .req {
  font-family: var(--r-code-font);
  font-size: 0.8em;
  font-weight: 700;
  color: var(--accent);
  background: var(--req-bg);
  border: 1px solid var(--line);
  border-radius: 5px;
  padding: 0.02em 0.35em;
  margin-right: 0.3em;
}

/* KPI row for headline numbers */
.reveal.deck-theme .kpis {
  display: grid;
  grid-template-columns: repeat(var(--n, 3), 1fr);
  gap: 0.7em;
  margin-top: 0.4em;
}

/* A two-column "figure + commentary" layout */
.reveal.deck-theme .fig-split {
  display: grid;
  grid-template-columns: var(--cols, 1.15fr 1fr);
  gap: 1em;
  align-items: center;
}
.reveal.deck-theme .fig-split .figure { width: 100%; text-align: center; }
/* .figure img/video are display:block with margin:0, so a picture narrower than
   its (100%-wide) split column would hug the left edge — text-align can't move
   a block. Recentre with auto side margins; a picture that fills the column is
   unaffected. */
.reveal.deck-theme .fig-split .figure img,
.reveal.deck-theme .fig-split .figure video { margin-inline: auto; }
.reveal.deck-theme .fig-split .figure-placeholder { width: 100%; }

/* Rounded portrait chips, the same trick DeckBrand uses for the logo bar: the
   corners come from a wrapper that clips its contents, not from the picture, so
   a headshot supplied as a flat card (name and role baked into the PNG, its own
   background colour and all) rounds off exactly like a bare portrait does.
   Sizing stays on the img — the chip shrink-wraps whatever it is given. */
.reveal.deck-theme .photo-chip {
  display: inline-flex;
  border-radius: 8px;
  overflow: hidden;
}
.reveal.deck-theme .photo-chip img { display: block; }

/* A centred, wrapping row of those chips — one per organisation on a team slide. */
.reveal.deck-theme .people {
  display: flex;
  justify-content: center;
  flex-wrap: wrap;
  gap: 0.4em;
  margin-top: 0.3em;
}

/* Stand-in for a not-yet-sourced image: a dashed frame naming exactly what to
   shoot/screenshot and the filename it will slot into once dropped in place. */
.reveal.deck-theme .figure-placeholder {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 0.4em;
  background: var(--surface);
  border: 2px dashed var(--line);
  border-radius: 12px;
  padding: 0.8em 1.1em;
  min-height: 170px;
  text-align: center;
}
.reveal.deck-theme .figure-placeholder.compact {
  min-height: 92px;
  padding: 0.5em 1em;
}
.reveal.deck-theme .figure-placeholder__tag {
  font-family: var(--r-code-font);
  text-transform: uppercase;
  letter-spacing: 0.12em;
  font-size: 0.42em;
  color: var(--accent-orange);
  border: 1px solid var(--accent-orange);
  border-radius: 999px;
  padding: 0.15em 0.9em;
}
.reveal.deck-theme .figure-placeholder__desc {
  font-size: 0.56em;
  color: var(--r-main-color);
  max-width: 34em;
  line-height: 1.3;
}
.reveal.deck-theme .figure-placeholder__target {
  font-family: var(--r-code-font);
  font-size: 0.44em;
  color: var(--comment);
}
</style>
