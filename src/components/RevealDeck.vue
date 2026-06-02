<script setup>
import { ref, onMounted, onBeforeUnmount, nextTick } from 'vue'
import Reveal from 'reveal.js'
import 'reveal.js/reveal.css'

// A thin Vue wrapper around reveal.js. Slot content must be a sequence of
// <section> elements (the slides). The deck is initialised on mount and torn
// down on unmount so navigating away leaves no global key/resize listeners.
const props = defineProps({
  options: { type: Object, default: () => ({}) }
})

const root = ref(null)
let deck = null

onMounted(async () => {
  await nextTick()
  deck = new Reveal(root.value, {
    // We run inside a vue-router hash route, so reveal's own hash navigation
    // would fight with the router — keep it off and drive with the keyboard.
    hash: false,
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
.reveal.deck-theme .center { text-align: center; }

/* Title slide */
.reveal.deck-theme .title-slide h1 {
  font-size: 2.4em;
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
  margin-top: 1.2em;
  color: var(--comment);
  font-size: 0.7em;
}
.reveal.deck-theme .title-slide .venue-note {
  display: inline-block;
  margin-top: 1.4em;
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
.reveal.deck-theme .figure img {
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
</style>
