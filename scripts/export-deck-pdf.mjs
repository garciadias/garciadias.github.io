#!/usr/bin/env node
/**
 * Export a reveal.js deck from the built site to a print-quality PDF.
 *
 * Renders reveal's own print view (?print-pdf) so each page is exactly the
 * slide as it appears on screen, forces the light theme, strips the app's
 * on-screen chrome, and scales each page up to 3840x2160 px (4K) with a CSS
 * transform, so nothing re-flows and the text/SVG stay vector — sharp at any
 * zoom, with the page box at a true 4K size.
 *
 * Usage:
 *   node scripts/export-deck-pdf.mjs <deck-id> [outfile] [--width=3840]
 *                                    [--url=http://localhost:4180]
 *                                    [--no-slide-numbers]
 *
 * Requires a static server for dist/ to already be running (`npx vite preview`).
 */
import { execFile } from 'node:child_process'
import { existsSync } from 'node:fs'
import { promisify } from 'node:util'
import puppeteer from 'puppeteer-core'

const run = promisify(execFile)

const argv = process.argv.slice(2)
const flag = (name, fallback) => {
  const hit = argv.find(a => a.startsWith(`--${name}=`))
  return hit ? hit.split('=').slice(1).join('=') : fallback
}
const has = name => argv.includes(`--${name}`)
const positional = argv.filter(a => !a.startsWith('--'))

const deckId = positional[0]
if (!deckId) {
  console.error('usage: node scripts/export-deck-pdf.mjs <deck-id> [outfile] [--width=3840]')
  process.exit(1)
}
const outFile = positional[1] || `${deckId}.pdf`
const origin = flag('url', 'http://localhost:4180').replace(/\/$/, '')
const targetWidth = Number(flag('width', 3840))
const targetHeight = Math.round(targetWidth * 9 / 16)

const chrome = ['/usr/bin/google-chrome', '/usr/bin/google-chrome-stable', '/usr/bin/chromium', '/snap/bin/chromium']
  .find(p => existsSync(p))
if (!chrome) throw new Error('no Chrome/Chromium binary found')

// ?print-pdf must sit in the query string (reveal reads location.search); the
// deck route lives in the hash because vue-router uses hash history here.
const url = `${origin}/?print-pdf#/presentations/${deckId}`

const browser = await puppeteer.launch({
  executablePath: chrome,
  headless: true,
  args: ['--no-sandbox', '--disable-gpu', '--font-render-hinting=none', '--force-device-scale-factor=1']
})

try {
  const page = await browser.newPage()

  // Light theme: the app reads localStorage first, then prefers-color-scheme.
  // Pin both so the export never inherits the machine's dark preference.
  await page.emulateMediaFeatures([{ name: 'prefers-color-scheme', value: 'light' }])
  await page.evaluateOnNewDocument(() => {
    try { localStorage.setItem('theme', 'light') } catch { /* private mode */ }
  })

  await page.setViewport({ width: 1600, height: 900, deviceScaleFactor: 1 })
  await page.goto(url, { waitUntil: 'networkidle0', timeout: 90_000 })
  await page.waitForSelector('.pdf-page', { timeout: 60_000 })
  await page.evaluate(() => document.fonts.ready)
  await new Promise(r => setTimeout(r, 2000))

  // Reveal writes the print page size onto body as an inline style, derived
  // from its width/height config plus the configured margin.
  const slide = await page.evaluate(() => ({
    w: parseFloat(document.body.style.width),
    h: parseFloat(document.body.style.height)
  }))
  if (!slide.w || !slide.h) throw new Error(`could not read reveal print page size (${slide.w}x${slide.h})`)

  // Blow the slide box up to the target page size. A transform (rather than a
  // zoom or a font-size change) means nothing re-flows: the page is the exact
  // layout reveal computed, drawn larger, and still vector.
  const scale = Math.min(targetWidth / slide.w, targetHeight / slide.h)

  // The deck view is a `position: fixed` fullscreen overlay. Fixed elements
  // don't paginate, so reveal's print pages would collapse into a single
  // printed page. Lift .reveal out to be body's only child and let them flow.
  const pages = await page.evaluate(({ w, h, pageW, pageH, scale, hideNumbers }) => {
    const reveal = document.querySelector('.reveal')
    reveal.remove()
    document.body.replaceChildren(reveal)
    document.body.style.cssText = `margin:0;padding:0;width:${pageW}px;height:auto;overflow:visible;background:#fff`
    document.documentElement.style.cssText = 'height:auto;overflow:visible'

    const style = document.createElement('style')
    style.textContent = `
      /* Reveal injects its own @page sized to the slide box; ours must win or
         each scaled-up page would be split across several sheets. */
      @page { size: ${pageW}px ${pageH}px; margin: 0; }
      /* The app's overlay controls are screen-only furniture. */
      .deck-chrome, .laser-dot, .reveal .controls, .reveal .progress { display: none !important; }
      ${hideNumbers ? '.slide-number, .slide-number-pdf { display: none !important; }' : ''}
      /* Let the print pages stack and break instead of being clipped. */
      .reveal, .reveal .slides {
        position: static !important; width: auto !important; height: auto !important;
        overflow: visible !important; transform: none !important; left: auto !important; top: auto !important;
      }
      .reveal .pdf-page {
        position: relative !important; width: ${pageW}px !important; height: ${pageH}px !important;
        page-break-after: always; break-after: page; overflow: hidden;
      }
      .reveal .pdf-page:last-child { page-break-after: auto; break-after: auto; }
    `
    document.head.appendChild(style)

    const pdfPages = [...document.querySelectorAll('.pdf-page')]
    for (const p of pdfPages) {
      const scaler = document.createElement('div')
      scaler.style.cssText =
        `position:absolute;left:0;top:0;width:${w}px;height:${h}px;` +
        `transform:scale(${scale});transform-origin:0 0`
      scaler.append(...p.childNodes)
      p.appendChild(scaler)
    }
    return pdfPages.length
  }, {
    w: slide.w, h: slide.h,
    pageW: targetWidth, pageH: targetHeight,
    scale, hideNumbers: has('no-slide-numbers')
  })

  if (!pages) throw new Error('reveal produced no .pdf-page elements')
  await new Promise(r => setTimeout(r, 1500))

  await page.pdf({
    path: outFile,
    width: `${targetWidth}px`,
    height: `${targetHeight}px`,
    margin: { top: 0, right: 0, bottom: 0, left: 0 },
    printBackground: true,
    preferCSSPageSize: false,
    pageRanges: `1-${pages}`
  })

  const { stdout } = await run('pdfinfo', [outFile]).catch(() => ({ stdout: '' }))
  const size = stdout.match(/Page size:\s*(.+)/)?.[1]?.trim() ?? ''
  console.log(
    `${outFile} — ${pages} pages, ${size} (= ${targetWidth}x${targetHeight} px @96dpi), ` +
    `vector, light theme, slide box ${slide.w}x${slide.h} scaled ${scale.toFixed(3)}x`
  )
} finally {
  await browser.close()
}
