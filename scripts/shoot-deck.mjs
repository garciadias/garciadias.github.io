#!/usr/bin/env node
/**
 * Screenshot every slide of a deck plus a per-slide overflow report.
 *
 * Usage: node scripts/shoot-deck.mjs <deck-id> <outdir> [--url=http://localhost:5173]
 */
import { mkdirSync, writeFileSync } from 'node:fs'
import puppeteer from 'puppeteer-core'

const argv = process.argv.slice(2)
const flag = (n, d) => {
  const hit = argv.find(a => a.startsWith(`--${n}=`))
  return hit ? hit.split('=').slice(1).join('=') : d
}
const pos = argv.filter(a => !a.startsWith('--'))
const deckId = pos[0]
const outDir = pos[1]
const origin = flag('url', 'http://localhost:5173').replace(/\/$/, '')
if (!deckId || !outDir) {
  console.error('usage: node scripts/shoot-deck.mjs <deck-id> <outdir> [--url=...]')
  process.exit(1)
}
mkdirSync(outDir, { recursive: true })

const browser = await puppeteer.launch({
  executablePath: process.env.CHROME_BIN || '/snap/bin/chromium',
  headless: true,
  args: ['--no-sandbox', '--disable-gpu', '--force-device-scale-factor=1', '--hide-scrollbars']
})
const page = await browser.newPage()
await page.setViewport({ width: 1280, height: 720, deviceScaleFactor: 1.5 })

// Pin the theme, defaulting to light — that is what the deck presents in (see
// the LIGHT block in RevealDeck.vue and scripts/export-deck-pdf.mjs), and the
// app reads localStorage first, then prefers-color-scheme, so both need
// setting or the shot inherits the machine's dark preference.
const theme = flag('theme', 'light')
await page.emulateMediaFeatures([{ name: 'prefers-color-scheme', value: theme }])
await page.evaluateOnNewDocument((t) => {
  try { localStorage.setItem('theme', t) } catch { /* private mode */ }
}, theme)

// Present mode (?deck=<id>): reveal owns the hash there, so slides can be
// addressed directly as #/<n> without fighting vue-router.
await page.goto(`${origin}/?deck=${deckId}`, { waitUntil: 'networkidle2', timeout: 90000 })
await page.waitForSelector('.reveal .slides > section', { timeout: 60000 })
await new Promise(r => setTimeout(r, 3000))

const count = await page.evaluate(() => document.querySelectorAll('.reveal .slides > section').length)
const report = []

for (let i = 0; i < count; i++) {
  await page.evaluate(n => { window.location.hash = `/${n}` }, i)
  await new Promise(r => setTimeout(r, 900))
  const info = await page.evaluate(() => {
    const s = document.querySelector('.reveal .slides > section.present') ||
      document.querySelector('.reveal .slides > section')
    const wrap = document.querySelector('.reveal .slides')
    const m = (wrap?.style.transform || '').match(/scale\(([\d.]+)\)/)
    const scale = m ? Number(m[1]) : 1
    const wrapH = wrap ? wrap.clientHeight : 0
    const r = s.getBoundingClientRect()
    const imgs = [...s.querySelectorAll('img')].map(im => ({
      src: (im.getAttribute('src') || '').split('/').pop(),
      ok: im.complete && im.naturalWidth > 0
    }))
    // Deepest element that pokes outside the slide box, if any.
    const sb = s.getBoundingClientRect()
    let worst = null
    for (const el of s.querySelectorAll('*')) {
      const b = el.getBoundingClientRect()
      if (!b.width && !b.height) continue
      const over = Math.max(b.bottom - sb.bottom, sb.top - b.top, b.right - sb.right, sb.left - b.left)
      if (over > 4 && (!worst || over > worst.over)) {
        worst = { over: Math.round(over), tag: el.tagName.toLowerCase(), cls: el.className?.toString?.().slice(0, 60) || '' }
      }
    }
    // Slide content that lands on top of the persistent DeckBrand chips. The
    // brand bar is anchored to the .reveal box, not the slide, so a tall slide
    // simply runs over it; .slide-body exists to veil text that does. Report
    // the deepest text node that overlaps, so collisions are measured rather
    // than eyeballed.
    const chips = [...document.querySelectorAll('.deck-brand-bottom .brand-chip, .deck-brand-top-right')]
      .map(c => c.getBoundingClientRect())
      .filter(b => b.width > 0 && b.height > 0)
    // Measure where the glyphs actually are, not the element box: a centred
    // full-width <p> overlaps the bottom-left chip on paper while its text sits
    // harmlessly in the middle of the slide. Range rects give one box per
    // rendered line, so centred and short lines report their true extent.
    const inkRects = (el) => {
      const out = []
      for (const node of el.childNodes) {
        if (node.nodeType !== Node.TEXT_NODE || !node.textContent.trim()) continue
        const r = document.createRange()
        r.selectNodeContents(node)
        out.push(...r.getClientRects())
      }
      // Inline children (<strong>, <code>, <a>) carry ink too.
      for (const child of el.children) {
        if (getComputedStyle(child).display.startsWith('inline')) {
          out.push(...child.getClientRects())
        }
      }
      return out
    }
    const hits = []
    for (const el of s.querySelectorAll('p, li, h1, h2, h3, h4, td, th, .pill, .caption, .panel, img')) {
      if (el.tagName !== 'IMG' && !el.textContent.trim()) continue
      // Only leaf-ish text holders, so a container doesn't shadow its children.
      if (el.querySelector('p, li, h1, h2, h3, h4, .pill, .caption')) continue
      const boxes = el.tagName === 'IMG' || el.classList.contains('panel')
        ? [el.getBoundingClientRect()]
        : inkRects(el)
      let area = 0
      for (const b of boxes) {
        if (!b.width || !b.height) continue
        for (const c of chips) {
          const ox = Math.min(b.right, c.right) - Math.max(b.left, c.left)
          const oy = Math.min(b.bottom, c.bottom) - Math.max(b.top, c.top)
          if (ox > 2 && oy > 2) area += ox * oy
        }
      }
      if (area > 40) {
        hits.push({
          tag: el.tagName.toLowerCase(),
          area: Math.round(area),
          veiled: !!el.closest('.slide-body, .panel, .figure'),
          text: (el.tagName === 'IMG' ? el.getAttribute('src').split('/').pop()
            : el.textContent.trim().replace(/\s+/g, ' ')).slice(0, 40)
        })
      }
    }
    hits.sort((a, b) => b.area - a.area)

    const h = s.querySelector('h1,h2,h3')
    return {
      brandHits: hits.slice(0, 3),
      title: (h ? h.textContent : '').trim().replace(/\s+/g, ' ').slice(0, 70),
      scale: Number(scale.toFixed(3)),
      slideH: Math.round(s.scrollHeight),
      slideW: Math.round(s.scrollWidth),
      renderedH: Math.round(r.height),
      wrapH,
      overflowsViewport: Math.round(r.height - wrapH),
      broken: imgs.filter(x => !x.ok).map(x => x.src),
      worst,
      chars: s.innerText.replace(/\s+/g, ' ').trim().length
    }
  })
  const file = `${outDir}/slide-${String(i).padStart(2, '0')}.png`
  await page.screenshot({ path: file, clip: { x: 0, y: 0, width: 1280, height: 720 } })
  report.push({ i, file, ...info })
  const naked = info.brandHits.filter(x => !x.veiled)
  const bad = info.overflowsViewport > 6 || info.broken.length || info.worst || naked.length
  console.log(
    `[${String(i).padStart(2)}] ${bad ? 'CHECK' : 'ok   '} h=${info.renderedH}/${info.wrapH} ` +
    `scale=${info.scale} chars=${info.chars}` +
    (info.worst ? ` worst=+${info.worst.over}px <${info.worst.tag} ${info.worst.cls}>` : '') +
    (info.broken.length ? ` BROKEN:${info.broken.join(',')}` : '') +
    `  ${info.title}`
  )
  for (const x of info.brandHits) {
    console.log(`        ${x.veiled ? 'brand-overlap(veiled)' : 'BRAND COLLISION     '} ` +
      `${x.area}px² <${x.tag}> "${x.text}"`)
  }
}

writeFileSync(`${outDir}/report.json`, JSON.stringify(report, null, 2))
console.log(`\n${count} slides -> ${outDir}`)
await browser.close()
