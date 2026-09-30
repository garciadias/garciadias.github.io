<script setup>
import { onBeforeUnmount, onMounted, ref } from "vue";
import DeckBrand from "@/components/DeckBrand.vue";
import RevealDeck from "@/components/RevealDeck.vue";

// IAA-SO School on AI/ML in Astronomy 2026 · Unsupervised Learning pillar.
// A ~90-minute lecture. Narrative: stars born together share a chemical
// fingerprint, and clusters dissolve over time — so Galactic archaeology
// wants to *reconstruct* them from abundances alone (chemical tagging,
// Freeman & Bland-Hawthorn 2002; Kos et al. 2017). That is a clustering
// problem in a high-D, noisy space, so the talk is a tour of the clustering
// toolbox, each algorithm chosen to fix the failure mode of the last:
//
//   K-means (centroid, needs K) → KNN (the nearest-neighbour primitive) →
//   DBSCAN (density,
//   arbitrary shape) → HDBSCAN* (density at all scales) → PLSCAN
//   (persistence, no min-cluster-size) → t-SNE (project then see) → UMAP
//   (graph embedding, global) → EVoC (fuses UMAP + HDBSCAN* + PLSCAN).
//
// Source material: the Garcia-Dias et al. 2020 "Clustering analysis" chapter
// (K-means, GMM, DBSCAN, Ward; SSE/silhouette/dip/homogeneity — note it is a
// psychiatry/neuroimaging methods chapter, not an astronomy paper, and it does
// not cover KNN), Bot et al. 2025
// (PLSCAN, persistence-based multiscale density clustering), and Kos et al.
// 2017 (t-SNE chemical tagging) — then the benchmark we run in this repo:
// t-SNE vs UMAP vs EVoC on APOGEE DR19 + Gaia DR3, scored against kinematic
// ground truth.
//
// Style follows the PyData lightning deck (same RevealDeck + DeckBrand
// skeleton, same panel/checklist/crosslist vocabulary); the argument arc
// echoes the BigDataLDN talk ("why do some ML models fail" → no free lunch
// → know your data → never quit thinking), now pointed at stars.
// The speaker portrait and future career photos live in the shared deck assets
// (public/presentations/shared/, the same files the FLIP decks use), so the
// About-me slide needs both roots — same helper name as FlareDay2026Deck.
const asset = (name) =>
    `${import.meta.env.BASE_URL}presentations/iaa-so-chemical-tagging-2026/${name}`;
const sharedAsset = (name) =>
    `${import.meta.env.BASE_URL}presentations/${name}`;

// Brand marks live in the site's shared logo store — the same directory
// DeckBrand reads the QR chip from — rather than in the deck's own asset
// folder, so the current-work strip is one line per logo and stays in step
// with the site if a mark is ever replaced there.
const logoAsset = (name) =>
    `${import.meta.env.BASE_URL}static/img/logos/${name}`;

// ── The full-bleed video slide (0c) ──────────────────────────────────────────
// A star-cluster formation simulation, edge to edge, started the moment the
// slide is reached so the audience sees the collapse in about 45 seconds.
// reveal's own data-background-* hooks cannot be asked for a playback rate —
// they hand the URL straight to an iframe — so the video gets its own player,
// built by YouTube's IFrame API (which also handles the postMessage handshake a
// raw `enablejsapi=1` embed never performs: with the handshake missing, the
// player is silent about being ready and no command can be delivered at all).
// Playback, rate and the unmute happen when the section gains reveal's
// `present` class and pause when it loses it, so audio never leaks onto the
// next slide. The rate asked for below is 3x because that is what the talk
// wanted; YouTube's player will not go past 2x, so 2x is what runs.
//
// The section is given the 1280x720 canvas size inline on purpose: reveal's
// `center: true` sets `min-height: 0`, so an empty section collapses to nothing
// and the player would have no box to fill.
const VIDEO = {
    id: "3z9ZKAkbMhY",
    // The talk wanted 3x; YouTube publishes the rates it honours and this video's
    // player offers 0.25 … 2, so a bare setPlaybackRate(3) is clamped back to 2
    // (measured: availablePlaybackRates = [0.25 … 2], playbackRate = 2 after the
    // call). The deck therefore asks for the fastest rate the player will deliver
    // and reports the gap once, rather than believing it got what it asked for.
    rate: 3,
    vars: {
        controls: 1,
        rel: 0,
        modestbranding: 1,
        playsinline: 1,
        iv_load_policy: 3,
        mute: 1,
    },
};

const videoSlide = ref(null);
const videoHost = ref(null);

let ytApi = null;
let player = null;
let playerReady = false;
let videoWanted = false;
let rateChecked = false;

const loadYouTubeApi = () =>
    (ytApi ||= new Promise((resolve) => {
        if (window.YT?.Player) return resolve(window.YT);
        // The API calls this global when it finishes loading; chain rather than
        // clobber, in case another deck on the page got there first.
        const prior = window.onYouTubeIframeAPIReady;
        window.onYouTubeIframeAPIReady = () => {
            prior?.();
            resolve(window.YT);
        };
        const script = document.createElement("script");
        script.src = "https://www.youtube.com/iframe_api";
        script.async = true;
        document.head.appendChild(script);
    }));

// YouTube publishes the rates it will honour for a given video and clamps
// anything outside that list, so ask for the fastest rate it can deliver.
const pickRate = () => {
    let rates = [];
    try {
        rates = player?.getAvailablePlaybackRates?.() || [];
    } catch {
        rates = [];
    }
    return rates.length
        ? (rates.filter((rate) => rate <= VIDEO.rate).pop() ?? rates[0])
        : VIDEO.rate;
};

// Unmuting a video the browser has not authorised to autoplay with sound stops
// it dead (measured on a link opened straight at this slide: the player reports
// buffering 3, then falls back to -1, the clock never leaves 0 and the big play
// button returns). So sound is asked for only when the page has actually been
// used — the key-press that advanced the deck counts — and the video otherwise
// plays silent and stays silent: promoting it to sound mid-flight pauses it, and
// the player then refuses a programmatic restart (also measured), so leaving it
// alone is the only predictable option. Stepping off the slide and back on is
// what brings the sound, since that path plays with activation behind it.
const soundAllowed = () => {
    const activation = navigator.userActivation;
    return activation ? activation.hasBeenActive : true;
};

const startVideo = () => {
    if (!playerReady) return;
    player.playVideo();
    player.setPlaybackRate(pickRate());
    if (soundAllowed()) player.unMute();
};

const stopVideo = () => {
    if (playerReady) player.pauseVideo();
};

const onPlayerStateChange = (event) => {
    // 1 = playing. Report the rate the player settled on, once: it should be the
    // one pickRate() asked for, and on this video that is deliberately not 3x.
    if (event.data !== 1 || rateChecked) return;
    rateChecked = true;
    const effective = player?.getPlaybackRate?.();
    const picked = pickRate();
    if (effective && Math.abs(effective - picked) > 0.01) {
        console.warn(
            `[iaa-so-chemical-tagging] video slide: asked for ${picked}x, player is running at ${effective}x`,
        );
    } else if (picked !== VIDEO.rate) {
        console.info(
            `[iaa-so-chemical-tagging] video slide: ${VIDEO.rate}x is what the notes ask for, ` +
                `${picked}x is the fastest this player offers`,
        );
    }
};

let videoObserver = null;
onMounted(async () => {
    const section = videoSlide.value;
    if (!section) return;

    // reveal drives the deck from the outside; the only signal it gives a slide
    // is its class list, so that is what we watch.
    const sync = () => {
        const present = section.classList.contains("present");
        if (present === videoWanted) return;
        videoWanted = present;
        // DeckBrand's chips are siblings of .slides, not children of the slide, so
        // no scoped rule can hide them for one slide — toggle them here instead.
        const chrome = document.querySelector(".deck-brand-bottom");
        if (chrome) chrome.style.display = present ? "none" : "";
        if (present) startVideo();
        else stopVideo();
    };
    videoObserver = new MutationObserver(sync);
    videoObserver.observe(section, {
        attributes: true,
        attributeFilter: ["class"],
    });
    sync();

    const YT = await loadYouTubeApi();
    // The deck can be left while the API script is still in flight.
    if (!videoHost.value) return;
    player = new YT.Player(videoHost.value, {
        videoId: VIDEO.id,
        host: "https://www.youtube-nocookie.com",
        playerVars: { ...VIDEO.vars, origin: window.location.origin },
        events: {
            onReady: () => {
                playerReady = true;
                // The API's replacement iframe carries width/height attributes, not
                // styles, and nothing in the deck CSS sizes it — fill the canvas here.
                const frame = player.getIframe?.();
                if (frame)
                    frame.style.cssText =
                        "display: block; width: 100%; height: 100%; border: 0";
                if (videoWanted) startVideo();
            },
            onStateChange: onPlayerStateChange,
            onError: (event) =>
                console.warn(
                    "[iaa-so-chemical-tagging] video slide: player error",
                    event?.data,
                ),
        },
    });
});

// ── The QR chip vs the slide surface ─────────────────────────────────────────
// The deck paints a near-opaque surface on .slides (see the style block at the
// end of this file), which would hide the QR chip: DeckBrand's chips are
// siblings of .slides, not children of any slide, so no scoped rule can reach
// them, and the theme keeps the chrome *under* the slides on purpose so that
// slide content can cover it. The chip is therefore lifted above the surface in
// CSS and steps aside on the slides whose content already occupies its corner —
// on those it was covered by that content before, so nothing changes visually,
// and every other slide keeps a QR the back of the room can still scan.
let chromeObserver = null;
let chromeRecheck = null;
let chromeSlides = null;
let chromeTransition = null;
let presentSlide = null;

const syncChrome = (section) => {
    const chrome = document.querySelector(".deck-brand-bottom");
    const chip = chrome?.querySelector(".deck-brand-group");
    if (!chrome || !chip || !section) return;
    // Measure the chip where it *would* sit. On every pass after the first the bar
    // is already display:none (either from this class or from the video slide), and
    // a hidden element measures as a 0×0 rect — which silently matches nothing and
    // would flip the class back off, showing the chip on top of the text it had
    // just stepped away from (measured: covered at 93ms, visible again at 793ms on
    // a jump from slide 0 to slide 60). The inline value is restored exactly, so
    // the video slide's own display:none is untouched.
    const inlineDisplay = chrome.style.display;
    chrome.style.display = "flex";
    const box = chip.getBoundingClientRect();
    chrome.style.display = inlineDisplay;
    const touches = (rect) =>
        rect &&
        rect.width > 1 &&
        rect.height > 1 &&
        rect.right > box.left + 2 &&
        rect.left < box.right - 2 &&
        rect.bottom > box.top + 2 &&
        rect.top < box.bottom - 2;
    // Ink, not boxes: a full-width svg or panel *box* reaches this corner on most
    // slides while its drawing does not, and box-testing hid the chip on 6 of 7
    // sampled slides. So measure what is actually drawn — the client rects of the
    // section's text nodes, plus any image (whose pixels cannot be inspected, so
    // it counts as ink everywhere inside its box, which is the honest reading).
    let covered = null;
    const walker = document.createTreeWalker(section, NodeFilter.SHOW_TEXT);
    for (
        let node = walker.nextNode();
        node && !covered;
        node = walker.nextNode()
    ) {
        if (!node.nodeValue?.trim()) continue;
        const range = document.createRange();
        range.selectNodeContents(node);
        if ([...range.getClientRects()].some(touches)) covered = "text";
        range.detach?.();
    }
    for (const element of section.querySelectorAll(
        "img, svg text, canvas, video",
    )) {
        if (touches(element.getBoundingClientRect()))
            covered = element.tagName.toLowerCase();
    }
    chrome.classList.toggle("qr-chip-covered", Boolean(covered));
};

onMounted(() => {
    const slides = document.querySelector(".reveal .slides");
    if (!slides) return;
    const measure = () => {
        const section = document.querySelector(
            ".reveal .slides > section.present",
        );
        if (section) syncChrome(section);
    };
    const check = () => {
        const section = document.querySelector(
            ".reveal .slides > section.present",
        );
        // Fragments change classes too; only a slide change needs re-measuring.
        if (!section || section === presentSlide) return;
        presentSlide = section;
        syncChrome(section);
        // reveal moves the incoming slide while its class list is already settled,
        // so that first measurement happens mid-transition. Repeat it when the
        // transform lands, and again on a timer for the cases with no transition to
        // wait for (a jump straight to a slide, or an interrupted one) — the deck's
        // own QA script measures at 900ms, so the chip has to have settled by then.
        clearTimeout(chromeRecheck);
        chromeRecheck = setTimeout(measure, 700);
    };
    chromeObserver = new MutationObserver(check);
    chromeObserver.observe(slides, {
        subtree: true,
        attributes: true,
        attributeFilter: ["class"],
    });
    chromeSlides = slides;
    chromeTransition = measure;
    slides.addEventListener("transitionend", measure);
    check();
});

onBeforeUnmount(() => {
    chromeObserver?.disconnect();
    chromeObserver = null;
    clearTimeout(chromeRecheck);
    chromeSlides?.removeEventListener("transitionend", chromeTransition);
    chromeSlides = null;
    chromeTransition = null;
    videoObserver?.disconnect();
    videoObserver = null;
    try {
        player?.destroy?.();
    } catch {
        // destroying a player whose iframe is already gone throws — harmless
    }
    player = null;
    playerReady = false;
});
</script>

<template>
    <RevealDeck :options="{ center: true }" theme-class="astro-theme">
        <template #chrome>
            <DeckBrand :qr="asset('presentation.png')" :show-logos="false" />
        </template>

        <!-- 0 · Title -->
        <section class="title-slide center">
            <div class="eyebrow">
                IAA-SO School on AI/ML in Astronomy 2026 · Unsupervised Learning
            </div>
            <h1>Chemical tagging: finding lost star clusters</h1>
            <p class="subtitle">
                Stars born together share a chemical fingerprint. Clusters
                dissolve, chemistry doesn't, so we reconstruct them from
                abundances alone, working our way from
                <strong>K-means to EVoC</strong>.
            </p>
            <p class="venue-note">
                Rafael Garcia-Dias · IAA-SO 2026 · 90 min introduction · 90 min
                problem · 120 min hands-on
            </p>
            <p class="small muted">
                K-means → KNN → DBSCAN → HDBSCAN* → PLSCAN → t-SNE → UMAP → EVoC
            </p>
            <aside class="notes">
                (~2 min) One-sentence hook: a cluster's chemistry is a fossil
                that survives long after the cluster itself has scattered into
                the field. Frame the day's three blocks up front, the next slide
                lays them out in one picture: 90 minutes of toolbox (one
                algorithm at a time, each fixing the last one's failure, landing
                on EVoC), 90 minutes of the real problem, and 120 minutes where
                they run it themselves. Then we start.
            </aside>
        </section>

        <!-- 0a · How the day runs — the agenda, right after the title so the
         three blocks (90 + 90 + 120, the speaker's own plan for the school day)
         are known before any science starts. The download slide is named here
         by title because students otherwise meet it six slides later, mid-
         lecture; the panel wording mirrors the three dividers (0a → 45b → 49a),
         so the day is announced, then re-announced at each seam. -->
        <section>
            <div class="eyebrow">How today runs · five hours, three blocks</div>
            <h2>One arc, three blocks, and your laptop does the last one</h2>
            <div class="cols" style="--n: 3; margin-top: 0.5em">
                <div class="panel">
                    <h3>90 min · Introduction</h3>
                    <p class="small">
                        What chemical tagging is, and the eight methods that go
                        after it (K-means → EVoC) each one fixing the last one's
                        failure.
                    </p>
                </div>
                <div class="panel">
                    <h3>90 min · The problem</h3>
                    <p class="small">
                        The same question inside a real paper: can any of these
                        methods recover clusters we already know? It ends with
                        your assignment.
                    </p>
                </div>
                <div class="panel">
                    <h3>120 min · Hands-on</h3>
                    <p class="small">
                        You, on DR19: check the download, run your own cluster,
                        and send the result back as a pull request.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.5em">
                Early on there is a slide called
                <em>“Download the data while I talk”</em>, start it and forget
                it: it runs in the background for the rest of the lecture.
            </p>
            <aside class="notes">
                (~1 min) Give the day's shape before the science, so nobody has
                to guess when we stop and when they work. Say the three blocks
                and their lengths out loud, and name the download slide here: it
                arrives six slides later, in the middle of the introduction, and
                anyone who starts it now has 2.2 GB less to wait for at the end.
                Also say how the day finishes; they send a pull request to the
                school repository, so the workshop leaves a record rather than a
                laptop full of files.
            </aside>
        </section>

        <!-- 0b · About me. Before the science, because two of the papers on the
         reference slide are the speaker's own — the 2018 K-means run over the
         APOGEE spectra and the 2019 experiment this lecture re-creates on DR19
         — so the audience should know who is claiming that. Facts are the
         site's own record (src/content/experience.js, publications.json:
         UFRGS open clusters, IAC PhD on 250 000+ spectra, KCL/Neurofind, Monzo,
         Floe, AI Centre/KCL; 18 Q1 articles, 5 000+ citations, h-index 19), and
         the timeline is the same career arc the Boehringer deck drew, restated
         in the astro palette and ending at the current role rather than at a
         speculative next one. Deliberately photo-forward: the four cards are
         the slide's content and cards 2–4 (plus, optionally, the portrait) are
         slots to fill with the speaker's own material before the talk, which is
         why they carry the same dashed "slot in" frame the rest of the repo
         uses for a figure that has not been sourced yet. -->
        <section>
            <div class="eyebrow">About me · who is teaching this</div>
            <h2>From astrophysics to healthcare industry</h2>
            <div
                class="cols"
                style="
                    --n: 4;
                    gap: 0.7em;
                    margin-top: 0.55em;
                    align-items: start;
                "
            >
                <!-- 1 · Portrait — real image, but itself replaceable -->
                <div>
                    <div
                        class="figure"
                        style="
                            display: block;
                            aspect-ratio: 4 / 3;
                            min-height: 0;
                            padding: 0.3em;
                        "
                    >
                        <img
                            :src="asset('k-means.png')"
                            alt="Portrait of Rafael Garcia-Dias"
                            style="
                                display: block;
                                width: 100%;
                                height: 100%;
                                object-fit: contain;
                            "
                        />
                    </div>
                    <p class="caption center">
                        K-Means on APOGEE Publication (2018)
                    </p>
                </div>
                <!-- 1 · Portrait — real image, but itself replaceable -->
                <div>
                    <div
                        class="figure"
                        style="
                            display: block;
                            aspect-ratio: 4 / 3;
                            min-height: 0;
                            padding: 0.3em;
                        "
                    >
                        <img
                            :src="asset('apogee.jpg')"
                            alt="Portrait of Rafael Garcia-Dias"
                            style="
                                display: block;
                                width: 100%;
                                height: 100%;
                                object-fit: contain;
                            "
                        />
                    </div>
                    <p class="caption center">
                        Unsupervised learning on APOGEE (IAC, 2015-2018)
                    </p>
                </div>
                <!-- 1 · Portrait — real image, but itself replaceable -->
                <div>
                    <div
                        class="figure"
                        style="
                            display: block;
                            aspect-ratio: 4 / 3;
                            min-height: 0;
                            padding: 0.3em;
                        "
                    >
                        <img
                            :src="asset('iac.jpg')"
                            alt="Portrait of Rafael Garcia-Dias"
                            style="
                                display: block;
                                width: 100%;
                                height: 100%;
                                object-fit: contain;
                            "
                        />
                    </div>
                    <p class="caption center">
                        PhD in ML applied to astrophysics (IAC, 2015-2018)
                    </p>
                </div>
                <!-- 1 · Portrait — real image, but itself replaceable -->
                <div>
                    <div
                        class="figure"
                        style="
                            display: block;
                            aspect-ratio: 4 / 3;
                            min-height: 0;
                            padding: 0.3em;
                        "
                    >
                        <img
                            :src="asset('codemotion.jpg')"
                            alt="Portrait of Rafael Garcia-Dias"
                            style="
                                display: block;
                                width: 100%;
                                height: 100%;
                                object-fit: contain;
                            "
                        />
                    </div>
                    <p class="caption center">Codemotion Milan (2019)</p>
                </div>
            </div>

            <!-- career timeline, adapted from the boehringer deck's svg. colour is the
           astro palette doing its usual job — blue for the astronomy years, pink
           for healthcare ai, amber for industry ml — with the green ring marking
           the current role. layout facts, all measured rather than eyed: the
           slide's content box is the full 1280 px canvas (reveal 6 does not pad
           sections for its margin option), so one viewbox unit renders at ~1.12
           px and the box edges are hard clip edges. hence the track runs 215 →
           1040: the first node's label has to clear deckbrand's qr chip, which
           owns the bottom-left corner (measured: the chip's right edge lands at
           159 px, the first label starts at 172 px), and the last node's label
           has to keep a right margin — at x=1056 it ended 2.5 px from the clip
           edge and one render cut the "flip" off. node spacing
           is then set from each label's rendered width (longest is "founding ml
           engineer" at 157 units) so no two labels can touch. label baselines are
           hand-set — years above the bar, role + org below — because a
           text/anchor mismatch shows up as a clipped label, which no amount of
           css can rescue. -->
            <svg
                viewbox="0 0 1140 88"
                xmlns="http://www.w3.org/2000/svg"
                style="
                    width: 100%;
                    height: auto;
                    display: block;
                    margin-top: 0.6em;
                "
                role="img"
                aria-label="career timeline. 2007–15 physics and open clusters at ufrgs, brazil. 2015–18 phd in machine learning on apogee spectra at the iac, tenerife. 2018–22 research associate and founding engineer at king's college london. 2022–23 decision scientist in credit risk at monzo. 2023–24 founding machine learning engineer at floe oral care. 2024 to now senior ai engineer at the ai centre and king's college london."
            >
                <!-- track bar, starting at the first node rather than at the slide edge -->
                <rect
                    x="215"
                    y="34"
                    width="825"
                    height="5"
                    rx="2.5"
                    style="fill: var(--line)"
                />

                <!-- 1 · ufrgs · physics, open clusters -->
                <circle
                    cx="215"
                    cy="36.5"
                    r="7.5"
                    style="fill: var(--accent-cyan)"
                    stroke="#ffffff"
                    stroke-width="2"
                />
                <line
                    x1="215"
                    y1="22"
                    x2="215"
                    y2="34"
                    style="stroke: var(--accent-cyan)"
                    stroke-width="1.4"
                />
                <text
                    x="215"
                    y="16"
                    text-anchor="middle"
                    style="
                        fill: var(--comment);
                        font-family: var(--r-code-font);
                        font-size: 12px;
                        font-weight: 700;
                    "
                >
                    2007–15
                </text>
                <text
                    x="215"
                    y="58"
                    text-anchor="middle"
                    style="
                        fill: var(--r-main-color);
                        font-size: 14.5px;
                        font-weight: 600;
                    "
                >
                    Physics · clusters
                </text>
                <text
                    x="215"
                    y="73"
                    text-anchor="middle"
                    style="fill: var(--comment); font-size: 13px"
                >
                    UFRGS · Brazil
                </text>

                <!-- 2 · IAC PhD — the APOGEE years -->
                <circle
                    cx="376"
                    cy="36.5"
                    r="7.5"
                    style="fill: var(--accent)"
                    stroke="#ffffff"
                    stroke-width="2"
                />
                <line
                    x1="376"
                    y1="22"
                    x2="376"
                    y2="34"
                    style="stroke: var(--accent)"
                    stroke-width="1.4"
                />
                <text
                    x="376"
                    y="16"
                    text-anchor="middle"
                    style="
                        fill: var(--comment);
                        font-family: var(--r-code-font);
                        font-size: 12px;
                        font-weight: 700;
                    "
                >
                    2015–18
                </text>
                <text
                    x="376"
                    y="58"
                    text-anchor="middle"
                    style="
                        fill: var(--r-main-color);
                        font-size: 14.5px;
                        font-weight: 600;
                    "
                >
                    PhD · ML on spectra
                </text>
                <text
                    x="376"
                    y="73"
                    text-anchor="middle"
                    style="fill: var(--comment); font-size: 13px"
                >
                    IAC · Tenerife
                </text>

                <!-- 3 · KCL — founding engineer on clinical ML -->
                <circle
                    cx="543"
                    cy="36.5"
                    r="7.5"
                    style="fill: var(--accent-pink)"
                    stroke="#ffffff"
                    stroke-width="2"
                />
                <line
                    x1="543"
                    y1="22"
                    x2="543"
                    y2="34"
                    style="stroke: var(--accent-pink)"
                    stroke-width="1.4"
                />
                <text
                    x="543"
                    y="16"
                    text-anchor="middle"
                    style="
                        fill: var(--comment);
                        font-family: var(--r-code-font);
                        font-size: 12px;
                        font-weight: 700;
                    "
                >
                    2018–22
                </text>
                <text
                    x="543"
                    y="58"
                    text-anchor="middle"
                    style="
                        fill: var(--r-main-color);
                        font-size: 14.5px;
                        font-weight: 600;
                    "
                >
                    Research associate
                </text>
                <text
                    x="543"
                    y="73"
                    text-anchor="middle"
                    style="fill: var(--comment); font-size: 13px"
                >
                    KCL · Neurofind
                </text>

                <!-- 4 · Monzo — credit-risk ML -->
                <circle
                    cx="702"
                    cy="36.5"
                    r="7.5"
                    style="fill: var(--accent-orange)"
                    stroke="#ffffff"
                    stroke-width="2"
                />
                <line
                    x1="702"
                    y1="22"
                    x2="702"
                    y2="34"
                    style="stroke: var(--accent-orange)"
                    stroke-width="1.4"
                />
                <text
                    x="702"
                    y="16"
                    text-anchor="middle"
                    style="
                        fill: var(--comment);
                        font-family: var(--r-code-font);
                        font-size: 12px;
                        font-weight: 700;
                    "
                >
                    2022–23
                </text>
                <text
                    x="702"
                    y="58"
                    text-anchor="middle"
                    style="
                        fill: var(--r-main-color);
                        font-size: 14.5px;
                        font-weight: 600;
                    "
                >
                    Decision scientist
                </text>
                <text
                    x="702"
                    y="73"
                    text-anchor="middle"
                    style="fill: var(--comment); font-size: 13px"
                >
                    Monzo · credit risk
                </text>

                <!-- 5 · Floe — founding ML engineer -->
                <circle
                    cx="870"
                    cy="36.5"
                    r="7.5"
                    style="fill: var(--accent-orange)"
                    stroke="#ffffff"
                    stroke-width="2"
                />
                <line
                    x1="870"
                    y1="22"
                    x2="870"
                    y2="34"
                    style="stroke: var(--accent-orange)"
                    stroke-width="1.4"
                />
                <text
                    x="870"
                    y="16"
                    text-anchor="middle"
                    style="
                        fill: var(--comment);
                        font-family: var(--r-code-font);
                        font-size: 12px;
                        font-weight: 700;
                    "
                >
                    2023–24
                </text>
                <text
                    x="870"
                    y="58"
                    text-anchor="middle"
                    style="
                        fill: var(--r-main-color);
                        font-size: 14.5px;
                        font-weight: 600;
                    "
                >
                    Founding ML engineer
                </text>
                <text
                    x="870"
                    y="73"
                    text-anchor="middle"
                    style="fill: var(--comment); font-size: 13px"
                >
                    Floe Oral Care
                </text>

                <!-- 6 · AI Centre / KCL — the current role, ringed -->
                <circle
                    cx="1040"
                    cy="36.5"
                    r="12"
                    fill="none"
                    style="stroke: var(--accent-green)"
                    stroke-width="1.6"
                />
                <circle
                    cx="1040"
                    cy="36.5"
                    r="7.5"
                    style="fill: var(--accent-pink)"
                    stroke="#ffffff"
                    stroke-width="2"
                />
                <line
                    x1="1040"
                    y1="22"
                    x2="1040"
                    y2="34"
                    style="stroke: var(--accent-pink)"
                    stroke-width="1.4"
                />
                <text
                    x="1040"
                    y="16"
                    text-anchor="middle"
                    style="
                        fill: var(--comment);
                        font-family: var(--r-code-font);
                        font-size: 12px;
                        font-weight: 700;
                    "
                >
                    2024→now
                </text>
                <text
                    x="1040"
                    y="58"
                    text-anchor="middle"
                    style="
                        fill: var(--r-main-color);
                        font-size: 14.5px;
                        font-weight: 600;
                    "
                >
                    Senior AI engineer
                </text>
                <text
                    x="1040"
                    y="73"
                    text-anchor="middle"
                    style="fill: var(--comment); font-size: 13px"
                >
                    AIC / KCL · FLIP
                </text>
            </svg>
            <aside class="notes">
                (~90 s) Keep this brisk; it earns its place on three points.
                One, why you are listening to me: two of the papers on tonight's
                reference slide are mine, and the 2019 one is the baseline we
                spend the back half of the talk re-creating on DR19, so the
                benchmark you are about to see is partly a self-check, and I
                will say so when we get there. Two, the arc in the timeline,
                left to right: physics and open clusters at UFRGS (MSc on
                Galactic open clusters, SMC clusters and Pismis 24 from
                VVV-ESO), then the PhD at the IAC, which is where the K-means
                run over 250 000 APOGEE spectra came from. Three, what happened
                after astronomy: founding engineer on clinical ML at KCL
                (Neurofind, and the harmonisation tool Neuroharmony), a
                six-month turn in fintech credit risk at Monzo, a founding ML
                engineer role on a diagnostic product at Floe, and since 2024
                foundation models plus federated learning at the AI Centre and
                KCL. Point at the last two nodes when you say it: the
                production-ML craft is why the numbers in this talk arrive with
                seed-stability and batch-effect controls attached, and why the
                workshop ships a Docker path. The MONAI work runs alongside
                rather than as a separate job, core contributor, and chair of
                its federated-learning working group with NVIDIA, which is the
                bridge between the two halves of this slide. The four photo
                cards carry real material now, the 2018 paper's first page, the
                APOGEE survey, the IAC and Codemotion Milan; to swap any of
                them, drop a replacement into
                public/presentations/iaa-so-chemical-tagging-2026/ and change
                the name inside its asset(...).
            </aside>
        </section>

        <!-- 0b2 · Current work — the timeline above ends on "Senior AI
         engineer · AIC / KCL / FLIP", so this slide answers the obvious next
         question in the order the user asked for (FLIP, NHS, MONAI), with the
         article screenshot he supplied and the four brand marks. Facts are the
         site's own record (src/content/experience.js: FLIP as the Federated
         Learning & Interoperability Platform, launched publicly Feb 2026,
         GSTT/KCH plus the Thai partner BDMS, MONAI federated-learning working
         group chair co-chaired with NVIDIA, contributions to MONAI Core and
         Generative; the cross-continental study is in src/content/presentations
         .js) and the paper's own first page, which supplies the title and the
         co-first-author marks. Its margin watermark reads
         "arXiv:submit/8129764 [cs.LG] 25 Sep 2026" — a submission, not yet an
         announced arXiv ID — so the caption prints no arXiv number; it becomes
         public on 30 Sep 2026. Logos come from public/static/img/logos/: aic and
         flip were already there, nhs_logo.webp is copied from the FLIP repo's
         own UI assets, and monai_logo.png is the project's colour wordmark from
         docs/images/MONAI-logo-color.png in Project-MONAI/MONAI (the repo had
         no MONAI mark; the GitHub org avatar is an icon-only variant).
         NOTE the h2 says NHS, but BDMS is a Thai hospital group, not a Trust —
         the two Trusts in the federation are GSTT and KCH. -->
        <section class="dense">
            <div class="eyebrow">Current work · what I do now</div>
            <h2>
                Federated learning across the NHS, without moving patient data
            </h2>
            <div
                class="cols"
                style="
                    --n: 2;
                    grid-template-columns: 0.82fr 1.18fr;
                    gap: 0.9em;
                    margin-top: 0.45em;
                    align-items: start;
                "
            >
                <div>
                    <div
                        class="figure"
                        style="
                            display: block;
                            aspect-ratio: 653 / 509;
                            min-height: 0;
                            padding: 0.3em;
                        "
                    >
                        <img
                            :src="asset('FLIP_article.png')"
                            alt="First page of the FLIP cross-continental federated learning preprint"
                            style="
                                display: block;
                                width: 100%;
                                height: 100%;
                                object-fit: contain;
                            "
                        />
                    </div>
                    <p class="caption center">
                        Co-first author ·
                        <em
                            >“Making Cross-Continental Federated Learning
                            Repeatable with FLIP”</em
                        >
                        preprint, out 30 Sep 2026
                    </p>
                </div>
                <div>
                    <div class="panel" style="margin-bottom: 0.5em">
                        <h3>FLIP · the platform</h3>
                        <p class="small">
                            Core contributor to the open-source
                            <strong
                                >Federated Learning Interoperability
                                Platform</strong
                            >: one model trained across NHS Trusts (GSTT, KCH)
                            and a Thai partner (BDMS)
                            <strong
                                >patient data never leaves the hospital</strong
                            >.
                        </p>
                    </div>
                    <div class="panel" style="margin-bottom: 0.5em">
                        <h3>NHS · the deployment</h3>
                        <p class="small">
                            Production AI infrastructure and distributed
                            training across three large NHS Trusts, millions of
                            patients with KCL, GSTT, NVIDIA, deepc, Flower Labs
                            and OneLondon.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>MONAI · the standard</h3>
                        <p class="small">
                            <strong
                                >Federated Learning Chair of the MONAI Working
                                Group</strong
                            >, co-chairing with NVIDIA; contributor to MONAI
                            Core and Generative.
                        </p>
                    </div>
                </div>
            </div>
            <div
                class="cols"
                style="
                    --n: 4;
                    gap: 1.1em;
                    margin-top: 0.6em;
                    align-items: center;
                    max-width: 34em;
                    margin-left: auto;
                    margin-right: auto;
                "
            >
                <div>
                    <img
                        :src="logoAsset('aic_logo.png')"
                        alt="London AI Centre for Value-Based Healthcare"
                        style="
                            display: block;
                            width: 100%;
                            height: 2.3em;
                            object-fit: contain;
                        "
                    />
                </div>
                <div>
                    <img
                        :src="logoAsset('nhs_logo.webp')"
                        alt="NHS"
                        style="
                            display: block;
                            width: 100%;
                            height: 2.3em;
                            object-fit: contain;
                        "
                    />
                </div>
                <div>
                    <img
                        :src="logoAsset('flip_logo.png')"
                        alt="FLIP, the Federated Learning Interoperability Platform"
                        style="
                            display: block;
                            width: 100%;
                            height: 2.3em;
                            object-fit: contain;
                        "
                    />
                </div>
                <div>
                    <img
                        :src="logoAsset('monai_logo.png')"
                        alt="Project MONAI"
                        style="
                            display: block;
                            width: 100%;
                            height: 2.3em;
                            object-fit: contain;
                        "
                    />
                </div>
            </div>
            <aside class="notes">
                (~2 min) The self-description ends on "Senior AI engineer · AIC
                / KCL / FLIP", so this is the question a school audience asks
                next: what is the current work? Three sentences carry it, the
                platform (FLIP: privacy-preserving training across trusts, live
                since February 2026, the data never moves), the scale
                (production AI across NHS Trusts with NVIDIA, deepc, Flower Labs
                and OneLondon), and the standards work (chair of MONAI's
                federated-learning working group, alongside NVIDIA). The
                preprint is also the bridge back to this talk's own subject: it
                is federated learning in healthcare, but the question in it is
                the one they are about to spend the day on; how do you validate
                a model against something real? Say one honest thing about
                timing: if you are presenting on 30 September the caption is
                literally "out today"; before that, say "submitted, out this
                week", because the margin watermark on the page is a submission
                number, not an arXiv ID yet. The four marks under it are the
                answer to "who is this person working with", in one line: AI
                Centre, NHS, FLIP, MONAI.
            </aside>
        </section>

        <!-- 0c · Video: large star-cluster formation (djxatlanta, YouTube). The
         slide is nothing but the video: no eyebrow, no heading, no branding —
         the deck's chrome is hidden for the duration by the script above, since
         the QR chip and logo bar are siblings of .slides and cannot be
         suppressed from inside a slide. It starts when reached and pauses when
         left, at the fastest rate YouTube's player offers (2x for this video —
         the 3x the talk asked for is not on its list); the URL is
         youtu.be/3z9ZKAkbMhY. Read the script comment before changing
         anything here. -->
        <section
            ref="videoSlide"
            class="video-slide"
            style="width: 1280px; height: 720px; padding: 0"
        >
            <!-- The host is empty: YouTube's IFrame API replaces this div with its own
           iframe, so it is the one styled to fill the canvas after onReady (the
           API's own width/height attributes default to 640x390). -->
            <div ref="videoHost" style="width: 100%; height: 100%"></div>
            <aside class="notes">
                (~45 s, no talking) A simulation of a massive cluster forming: a
                giant molecular cloud collapses, fragments, and lights up. Let
                it play. Two sentences total: one before ("watch what a cluster
                looks like while it is being born; every star you see here
                formed from the same gas") and one after ("that is why their
                chemistry is a fingerprint"). It does the work the next slide
                then names. Speed: the footage is 1:23 and the deck asks for the
                fastest rate YouTube will give it, which for this video is 2x;
                about 42 seconds. If that is still too long for the room, the
                player's own controls are live: drag the progress bar, or drop
                the rate to 1.5x on the settings gear (YouTube offers 0.25 to 2,
                nothing faster; 3x is not on the list). Operationally: it starts
                when the slide is reached and pauses when you leave it, so audio
                never runs onto the next slide. If you open a shared link
                straight to this slide it starts silent, browsers refuse sound
                until the page has been clicked or typed in, and stepping off
                the slide and back on gives it sound (or use the player's own
                speaker button). In the talk itself you will have advanced the
                deck by keyboard, so it starts with sound. If it does not start
                at all, click it once; the player's controls are live. Coming
                back to the slide resumes where it stopped rather than
                restarting; to replay it from the top, click the player's own
                progress bar. Credit if asked: simulation by Matthew Bate
                (University of Exeter / UK Astrophysical Fluids Facility), as
                uploaded by djxatlanta, "Large Star Cluster Formation [720p]",
                youtu.be/3z9ZKAkbMhY.
            </aside>
        </section>

        <!-- 1 · The science question -->
        <section>
            <div class="eyebrow">The science question</div>
            <h2>Stars are born together, then drift apart</h2>
            <div class="cols" style="--n: 3; margin-top: 0.3em">
                <div class="panel">
                    <h3>The premise</h3>
                    <p class="small">
                        A cluster forms from one well-mixed cloud, so its
                        members are chemically homogeneous to better than 0.1
                        dex; an abundance pattern is a
                        <strong>fingerprint of common origin</strong>.
                    </p>
                </div>
                <div class="panel">
                    <h3>The dream</h3>
                    <p class="small">
                        <em>Galactic archaeology</em>: use that fingerprint to
                        reconstruct clusters that dissolved billions of years
                        ago, and to recover members flung far from the cluster
                        centre.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Why it is a clustering problem</h3>
                    <p class="small">
                        Nothing tells us how many clusters there are, or which
                        star belongs to which. We only have a cloud of points
                        that should hide
                        <strong>unlabelled overdensities</strong>.
                    </p>
                </div>
            </div>
            <p class="small muted center" style="margin-top: 0.5em">
                <a
                    href="https://arxiv.org/abs/1709.00794"
                    target="_blank"
                    rel="noopener"
                    >Kos et al. 2017</a
                >, after
                <a
                    href="https://ui.adsabs.harvard.edu/abs/2002ARA%26A..40..487F/abstract"
                    target="_blank"
                    rel="noopener"
                    >Freeman &amp; Bland-Hawthorn 2002</a
                >.
            </p>
            <aside class="notes">
                (~2 min) Motivation first, no ML yet. Say the definition out
                loud: chemical tagging is using the abundances measured in
                stellar atmospheres to reconstruct chemically homogeneous
                clusters that have already dispersed. Then land the third panel,
                because it is the hinge of the whole lecture: we have no labels,
                no cluster count and no shapes; only points. That is the
                definition of a clustering problem, and everything from here on
                is about which clustering algorithm survives contact with this
                particular cloud.
            </aside>
        </section>

        <!-- 2 · The data -->
        <section>
            <div class="eyebrow">The data</div>
            <h2>A 16-dimensional chemical space (C-space)</h2>
            <div
                class="fig-split"
                style="--cols: 1fr 1fr; margin-top: 0.3em; align-items: start"
            >
                <div>
                    <p class="small">
                        One star → one point in <strong>C-space</strong>, the
                        space of its elemental abundances relative to iron.
                        <strong>APOGEE DR19</strong> supplies the chemistry;
                        <strong>Gaia DR3</strong> astrometry is held back as
                        ground truth.
                    </p>
                    <ul class="dotlist small">
                        <li>16 abundance dimensions per star</li>
                        <li>~183 000 stars after quality cuts (SNR ≥ 100)</li>
                        <li>Each element standardised before any distance</li>
                    </ul>
                    <p class="small muted" style="margin-top: 0.4em">
                        C, N, O, Na, Mg, Al, Si, S, K, Ca, Ti, V, Cr, Mn and Ni
                        as [X/Fe], plus [Fe/H] itself. Rescaling every element
                        to zero median and unit σ stops the one with the widest
                        scatter from deciding every distance.
                    </p>
                </div>
                <div class="figure">
                    <img
                        :src="asset('cspace_corner.png')"
                        alt="Two abundance dimensions: field stars grey, cluster members magenta in one tight clump"
                        style="width: 100%; height: auto"
                    />
                </div>
            </div>
            <aside class="notes">
                (~2 min) Frame the input; this is the same dataset they will use
                in the hands-on. Make two points. First, standardisation is not
                cosmetic: every algorithm from here on is a statement about
                distance, and an unstandardised element with ten times the
                scatter would silently own that distance. Second, the plot is a
                trap; it shows 2 of 16 dimensions, and the red clump is only
                obvious because I coloured it. In the other 14 dimensions, and
                without the colours, nobody can see it. That gap is why the rest
                of the lecture exists.
            </aside>
        </section>

        <!-- 3 · The clustering problem -->
        <section>
            <div class="eyebrow">The problem</div>
            <h2>Finding groups in 16-D is not obvious</h2>
            <p class="small">
                In 2-D you can see an overdensity. In 16-D the volume explodes,
                "look for the peak" stops working, and clustering analysis earns
                its keep, under its own rules (Garcia-Dias et al. 2020).
            </p>
            <div class="cols" style="--n: 3; margin-top: 0.45em">
                <div class="panel">
                    <h3>Multimodality first</h3>
                    <p class="small">
                        An algorithm <em>always</em> returns groups, even from a
                        uniform cloud. Test for a multi-peaked distribution
                        first (dip test).
                    </p>
                </div>
                <div class="panel">
                    <h3>Four tasks</h3>
                    <p class="small">
                        Feature selection → similarity metric → grouping
                        criterion → validation. The grouping criterion is the
                        core.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>The trap</h3>
                    <p class="small">
                        "Distinct groups will be found even when there is no
                        multipeak distribution"; the interpretation is on you.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.5em">
                There are no bad models, only models used outside their
                assumptions.
            </p>
            <aside class="notes">
                (~2 min) The chapter's core warning, and the BigDataLDN opening
                move. K-means (and most algorithms) will happily hand you
                clusters from pure noise. The honest question (for chemical
                tagging too) is not "did I get clusters?" but "are they real?"
                That's why validation gets its own section later. This slide
                also seeds the "four tasks": every algorithm that follows
                differs mainly in task 2 (metric) and task 3 (grouping
                criterion).
            </aside>
        </section>

        <!-- 3b · Setup — the one slide the students act on during the talk: the
             ~2.2 GB hands-on download starts now and runs while the algorithms
             are taught, so the session at the end of the class begins from a
             warm cache. Commands are the published repository's own
             (day_4_clustering/README.md Quick start, docs/student_activities.md
             §0); ./run.sh is the wrapper that file ships (mode 100755 in the
             repo) and it only saves typing the mount flags. Everything below is
             present in the students' clone: verified against the org repo tree.
             The fork-and-PR strip is the standard flow for a repo they have no
             write access to — the repo is public with allow_forking true, and
             .github/ carries no PR template, while day4-tests.yml does run on
             pull_request for day_4_clustering/** paths (their CI needs a
             maintainer's approval on a first-time contributor's PR). -->
        <section class="dense">
            <div class="eyebrow">Hands-on · start this now</div>
            <h2>Download the data while I talk</h2>
            <p class="small">
                You need a <strong>GitHub account</strong>,
                <strong>git</strong> and <strong>Docker</strong> (with
                <strong>Compose</strong>) nothing else. git and Docker are
                section B of the School Installation Guide, so they are already
                on your laptop; Python, <code>uv</code> and every package live
                inside the container.
            </p>
            <pre
                class="small"
                style="text-align: left; max-width: 46em; margin: 0.5em auto"
            ><code># fork it on GitHub first: github.com/iaa-so-training/iaa-advanced-neural-networks-2026
git clone https://github.com/&lt;your-username&gt;/iaa-advanced-neural-networks-2026.git
cd iaa-advanced-neural-networks-2026/day_4_clustering
./run.sh download --all    # DR19 catalogue + embeddings · ~2.2 GB · once</code></pre>
            <p class="small muted center" style="margin-top: 0.3em">
                Windows: <code>.\run.ps1 download --all</code> &nbsp;·&nbsp; the
                plain <code>docker run</code> form (no wrapper) is in
                <code>day_4_clustering/README.md</code>
            </p>
            <div
                class="cols"
                style="--n: 2; margin-top: 0.5em; align-items: stretch"
            >
                <div class="panel">
                    <h3>It keeps going without you</h3>
                    <p class="small">
                        Leave it running and keep listening; the download is
                        resumable and skips whatever is already on disk, so a
                        dropped connection costs nothing, just re-run it. When
                        it lands, check the setup:
                        <code>./run.sh run --fast</code> (~2 min, no GPU).
                    </p>
                </div>
                <div class="panel flip">
                    <h3>At the end of the class</h3>
                    <p class="small">
                        We work through the hands-on from this same folder:
                        <code>./run.sh lab</code> →
                        <code>http://localhost:8889</code>, then open
                        <code>notebooks/chemical_tagging.ipynb</code>. The
                        walkthrough is <code>docs/student_activities.md</code>,
                        and your cluster is the one assigned to you in
                        <code>docs/cluster_assignment.md</code>.
                    </p>
                </div>
            </div>
            <p class="small center" style="margin-top: 0.45em">
                <strong>Send it back as a pull request:</strong> branch (<code
                    >git checkout -b my-cluster</code
                >), commit, push (<code>git push -u origin my-cluster</code>),
                then <strong>Contribute → Open pull request</strong> on your
                fork, base <code>main</code>.
            </p>
            <aside class="notes">
                (~90 s) The only slide they have to act on, so leave it up long
                enough to read. The download is ~2.2 GB, the 1.17 GB SDSS-V DR19
                catalogue plus the ~1.0 GB embeddings bundle from Hugging Face,
                and it does not need them while it runs: that is the whole point
                of doing it now, so it finishes during the algorithms. Say the
                prerequisite line plainly: git and Docker-with-Compose are the
                only things on their laptops; everything else is in the image
                and nothing is installed on the host. Two reassurances worth
                giving: the download is resumable and skips what is already on
                disk, so flaky wifi is not fatal; and the ~2 min smoke run
                (`./run.sh run --fast`) is how they know they are ready before
                the session starts. If anyone has no Docker at all, the native
                path is `uv sync` with Python ≥ 3.13, same commands and flags,
                see `day_4_clustering/README.md`. The notebook is JupyterLab on
                port **8889**, deliberately not Jupyter's usual 8888, which is
                often already taken on a school laptop, and the deck ships
                `.ipynb` files now, so say the port and the file name out loud
                once: it saves twenty people asking. The fork is not
                bureaucracy, say it in one line: they can only push to their own
                copy of the repo, so the fork is where their branch lives, click
                Fork on the school repo page first, then clone YOUR fork, which
                is why the clone line carries their username. Two things to tell
                them before they go: branch before they start editing, and the
                fork's Contribute → Open pull request button already targets
                base `main`. Expect a CI click: GitHub holds the workflow run on
                each student's *first* PR until a maintainer approves it, so
                their checks will sit pending until one of us hits Approve and
                run. The next slide is the fork page itself: show it here,
                before they clone.
            </aside>
        </section>

        <!-- 3c · Fork it — a deliberate repeat of 45d4 ("Then fork: your own copy
         to push to"), placed directly after the setup slide because that is where
         students are told to fork before cloning, so the button is on screen
         while they are asked to press it. Same h2, same image and same caption
         as 45d4 so the two read as one page; only the eyebrow, the lead line and
         the notes are written for this position (here the fork is the immediate
         instruction rather than the end of a reading sequence). If either slide
         is ever edited, edit both — this pair is one picture shown twice on
         purpose. -->
        <section class="dense">
            <div class="eyebrow">Hands-on · fork it now</div>
            <h2>Then fork: your own copy to push to</h2>
            <p class="small">
                Before you clone: the fork is what gives you a copy of your own,
                your branch, your push, and later your pull request. One click,
                and no permission needed from anyone.
            </p>
            <div style="text-align: center; margin-top: 0.5em">
                <div class="figure" style="padding: 0.3em">
                    <img
                        :src="asset('fork.png')"
                        alt="The school repository's front page: Public badge, Fork button, and the licence named in the sidebar"
                        style="display: block; max-height: 440px; width: auto"
                    />
                </div>
            </div>
            <p class="caption center">
                github.com/iaa-so-training/iaa-advanced-neural-networks-2026
            </p>
            <aside class="notes">
                (~1 min) The same screenshot as the later "read the basics"
                slide, repeated here on purpose so the fork button is on screen
                while they are being asked to use it. Point at three things in
                order: the Public badge; they can read and fork it without
                asking; the Fork button; that is the whole action; and the
                licence in the sidebar, which they will meet again later in the
                day. Then the line that saves questions: fork first, then clone
                YOUR fork as the previous slide says, because a branch only
                exists somewhere you can push it. Two reassurances if the room
                looks worried: the fork count on screen (1) shows this is a
                two-minute job, and if their own first pull request later sits
                with a pending check, that is expected, GitHub waits for a
                maintainer to approve a first-time contributor's run, and one of
                us will click it.
            </aside>
        </section>

        <!-- 4 · K-means — how it works -->
        <section>
            <div class="eyebrow">Tool one · K-means</div>
            <h2>K-means: the default tool</h2>
            <div
                class="fig-split"
                style="--cols: 1fr 1fr; margin-top: 0.3em; align-items: start"
            >
                <div>
                    <ol class="contribs small">
                        <li>
                            <strong>Choose K</strong> and scatter K centres
                            μ<sub>1</sub>…μ<sub>K</sub> at random.
                        </li>
                        <li>
                            <strong>Assign</strong> each star to its nearest
                            centre (Euclidean).
                        </li>
                        <li>
                            <strong>Update</strong> each centre to the mean of
                            its points.
                        </li>
                        <li>
                            <strong>Repeat</strong> 2–3 until no star changes
                            cluster.
                        </li>
                    </ol>
                    <p class="small" style="margin-top: 0.4em">
                        Steps 2–3 descend a single objective, the within-cluster
                        scatter
                        <code
                            >SSE = Σ<sub>j</sub> Σ<sub>i∈S<sub>j</sub></sub>
                            ‖x<sub>i</sub> − μ<sub>j</sub>‖²</code
                        >.
                    </p>
                    <p class="small muted" style="margin-top: 0.3em">
                        Partitional · centroid-based · hard labels; every star
                        ends up in exactly one cluster, and
                        <strong>K is an input</strong>, never a result.
                    </p>
                </div>
                <div class="figure">
                    <img
                        :src="asset('kmeans.gif')"
                        alt="K-means iterating: points recoloured by nearest centroid, centroids moving each step"
                        style="width: 100%; height: auto"
                    />
                </div>
            </div>
            <aside class="notes">
                (~3 min) First algorithm of the day, so walk all four steps out
                loud while the GIF loops, point at the screen: recolour, move,
                recolour, move, stop. Say it is from the 1950s (Steinhaus 1956;
                Lloyd 1982; MacQueen 1967) and still the default first thing
                anyone reaches for. Then plant the two seeds that the next
                slides harvest: (a) K is an input, not an output, the algorithm
                cannot tell you how many clusters exist; (b) "nearest centre in
                Euclidean distance" means the boundaries between clusters are
                straight lines, so K-means can only carve the space into convex
                cells. Both come back as failure modes shortly.
            </aside>
        </section>

        <!-- 5 · Choose K and scatter K centres -->
        <section>
            <div class="eyebrow">Tool one · K-means · step 1 of 4</div>
            <h2>Choose K and scatter K centres</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('kmeans_step1_init.png')"
                    alt="All points unassigned, K star-shaped centres dropped at random"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                K is an input: you decide how many clusters to look for, then
                drop K centres; here on K randomly chosen stars.
            </p>
            <aside class="notes">
                (~30 s) The only decision the algorithm cannot make for you.
                Land it: K is an input, never a result. Note the three centres
                already carry a colour and a shape; nothing is assigned yet, but
                each centre has an identity, so the colours on the next slide
                read as "belongs to that centre" rather than as three arbitrary
                groups.
            </aside>
        </section>

        <!-- 6 · Assign each star to its nearest centre -->
        <section>
            <div class="eyebrow">Tool one · K-means · step 2 of 4</div>
            <h2>Assign each star to its nearest centre</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('kmeans_step2_assign.png')"
                    alt="Points coloured by their nearest centre; centres have not moved"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Each line ties one star to the centre that owns it, and where
                two bundles meet, the boundary is a straight line.
            </p>
            <aside class="notes">
                (~30 s) Trace one line with a finger: that star, that centre,
                nothing else. Then step back and let the bundles do the work.
                Two things to name. The seam where two fans meet is dead
                straight, because "nearest centre" is decided by a perpendicular
                bisector; that geometry is the failure mode we come back to. And
                these spokes are long, which is the point: add up their squared
                lengths and you have the SSE the algorithm is trying to shrink.
            </aside>
        </section>

        <!-- 7 · Move each centre to the mean of its stars -->
        <section>
            <div class="eyebrow">Tool one · K-means · step 3 of 4</div>
            <h2>Move each centre to the mean of its stars</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('kmeans_step3_update.png')"
                    alt="Centres jumped to the mean of their assigned points; points have not moved"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Hollow marks show where each centre was; the arrow runs to the
                mean of the stars it owns. Not one star moved.
            </p>
            <aside class="notes">
                (~30 s) The update step, and the arrows are the whole point:
                from a bad start the centres travel a long way in one pass.
                Emphasise 'mean'; that is where the name comes from, and why a
                single outlier drags a centre. Say explicitly that the stars did
                not move and did not change colour; only the centres did.
            </aside>
        </section>

        <!-- 8 · Repeat until nothing changes -->
        <section>
            <div class="eyebrow">Tool one · K-means · step 4 of 4</div>
            <h2>Repeat until nothing changes</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('kmeans_step4_converged.png')"
                    alt="Final stable partition: centres at cluster means, no point switching"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Eight passes later the spokes are as short as they can get, and
                their total squared length <em>is</em> the SSE.
            </p>
            <aside class="notes">
                (~30 s) Hold this next to step 2: there the spokes reach right
                across the field, here each centre sits at the heart of a short,
                tidy fan. That shortening is the algorithm working, because the
                total squared spoke length is exactly the SSE, the two moves are
                just the two ways of shortening it, and each can only make it
                smaller, so the loop must stop. It stops at a local minimum,
                which is the hook for the slide after next.
            </aside>
        </section>

        <!-- 9 · The loop, end to end -->
        <section class="dense">
            <div class="eyebrow">Tool one · K-means · the loop</div>
            <h2>Two moves, repeated until nothing changes</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 880px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('kmeans_storyboard.png')"
                    alt="Six panels: the raw data, three centres dropped at random, the assignment drawn as spokes, the centres moving to their means, a second assignment, and the converged partition"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 880px; margin-inline: auto"
            >
                (c)→(d) and (e)→(f) are the same two moves. That repetition
                <em>is</em> the algorithm.
            </p>
            <aside class="notes">
                (~30 s) The recap the four step slides cannot give, because you
                can never see two of them at once. Point at (c) and (e): same
                move. At (d) and (f): same move. Everything K-means does is
                those two, alternated, until the colours stop changing. Then set
                up the next slide by pointing back at (b): every one of these
                panels was decided by where those three centres happened to
                land.
            </aside>
        </section>

        <!-- 10 · K-means — initialization & SSE -->
        <section>
            <div class="eyebrow">Tool one · K-means</div>
            <h2>Run it again; you may get a worse answer</h2>
            <div
                class="fig-split"
                style="
                    --cols: 1.05fr 0.95fr;
                    margin-top: 0.3em;
                    align-items: start;
                "
            >
                <div class="figure">
                    <img
                        :src="asset('kmeans_init.png')"
                        alt="Two K-means runs from different initialisations on the same data: a good init finds the three-way split (SSE 331), a bad init merges two clusters (SSE 553)"
                        style="width: 100%; height: auto"
                    />
                </div>
                <div>
                    <ul class="dotlist small">
                        <li>
                            <strong>Initialisation matters.</strong> Same data,
                            two seeds: the good run finds the three-way split
                            (SSE = 331), the bad run merges two clusters (SSE =
                            553) a worse <em>local</em> minimum.
                        </li>
                        <li>
                            <strong>SSE is not convex.</strong> K-means only
                            reaches a local minimum: in the chapter, 4 of 5 runs
                            score ≈89% against known labels, one collapses to
                            ≈49%.
                        </li>
                        <li>
                            <strong>SSE</strong> (sum of squared error) discards
                            those poor runs, but it barely separates the good
                            ones, and it can never choose K.
                        </li>
                        <li>
                            <strong>Guardrails:</strong> K-means++ seeding
                            (Arthur &amp; Vassilvitskii 2007) and many restarts
                            &gt;250 for one cortical segmentation (Nanetti et
                            al. 2009).
                        </li>
                    </ul>
                </div>
            </div>
            <aside class="notes">
                (~3 min) Two ideas, and the audience must not conflate them.
                First, the figure: same data, two different seeds, the good run
                finds the three-way split, the bad run merges two clusters into
                one at a visibly higher SSE. That worse answer is a local
                minimum, not a bug (cluster labels still carry no meaning, never
                compare cluster "2" across runs). Second, the real hazard: SSE
                is a non-convex objective, so K-means descends into whichever
                local minimum its seed is nearest, and the chapter's example has
                one run in five landing at half the accuracy of the others. The
                cheap fix is restarts plus keeping the lowest SSE; K-means is
                fast enough that hundreds of restarts still beat one run of
                anything smarter. Close by flagging the gap this leaves: SSE
                cannot choose K, because it falls monotonically as K grows.
                Silhouette (Rousseeuw 1987) is the usual tool for that, and we
                come back to it in the validation section, with the same caveat
                as everything else here: only believe a clear peak.
            </aside>
        </section>

        <!-- 11 · K-means limits -->
        <section class="denser">
            <div class="eyebrow">Tool one · K-means limits</div>
            <h2>It assumes round, equal-sized clusters</h2>
            <div
                class="fig-split"
                style="
                    --cols: 1.1fr 1fr;
                    margin-top: 0.35em;
                    align-items: start;
                "
            >
                <div class="figure" style="margin: 0">
                    <img
                        :src="asset('kmeans_fail.png')"
                        alt="K-means splitting the two half-moons and cutting the big cluster in two"
                        style="width: 100%; height: auto; display: block"
                    />
                </div>
                <div>
                    <ul class="crosslist small">
                        <li>
                            <strong>Non-spherical shapes</strong>, a half-moon
                            is split in two.
                        </li>
                        <li>
                            <strong>Different scales</strong>, one wide Gaussian
                            eats the rest.
                        </li>
                        <li>
                            <strong>Unbalanced sizes</strong>, the small group
                            is absorbed.
                        </li>
                        <li>
                            <strong>Straight boundaries</strong>, nearest centre
                            cuts with lines.
                        </li>
                    </ul>
                    <p class="small muted" style="margin-top: 0.4em">
                        Give K-means the right K and it still fails here.
                        Assigning each star to the nearest centre
                        <em>in Euclidean distance</em> can only carve the space
                        with straight edges. These are characteristics, not bugs
                        but they are exactly the geometry 16-D chemical space
                        throws at us: uneven, anisotropic, overlapping groups.
                    </p>
                </div>
            </div>
            <aside class="notes">
                (~2 min) Make it concrete with the chapter's own example (Fig.
                13.4): give K-means the correct K and it still fails when the
                groups violate its assumptions, non-spherical shapes, different
                scales, unbalanced sizes, and the straight decision boundaries
                that come with Euclidean distance to a centre. Say the line out
                loud: these are characteristics, not bugs. That is the setup for
                the whole tour; every algorithm that follows exists to remove
                one of these assumptions.
            </aside>
        </section>

        <!-- 12 · KNN -->
        <section class="dense">
            <div class="eyebrow">Tool two · KNN</div>
            <h2>KNN: the nearest-neighbour primitive</h2>
            <p class="small">
                K-nearest neighbours asks one question,
                <em>which k points are closest to x?</em> It is not a clustering
                algorithm; it is the <strong>local primitive</strong> that
                DBSCAN, HDBSCAN*, UMAP and EVoC all call.
            </p>
            <div
                class="fig-split"
                style="--cols: 1.8fr 1fr; margin-top: 0.25em"
            >
                <div class="panel">
                    <h3>The algorithm</h3>
                    <ol class="contribs tight small">
                        <li>
                            <strong>Pick k</strong>, the neighbourhood size, and
                            the method's only knob.
                        </li>
                        <li>
                            <strong>Measure</strong> the distance d(x,
                            x<sub>i</sub>) to every other point.
                        </li>
                        <li>
                            <strong>Keep</strong> the k smallest → the k nearest
                            neighbours of x.
                        </li>
                        <li>
                            <strong>Use them</strong>, a majority vote
                            (classification), the k-th distance (density), or
                            the edges themselves (graph).
                        </li>
                    </ol>
                </div>
                <div class="figure" style="margin: 0">
                    <img
                        :src="asset('knn_core_distance.png')"
                        alt="Core distance: the k-th nearest-neighbour radius is small in dense regions and large in sparse ones"
                        style="width: 100%; height: auto; display: block"
                    />
                </div>
            </div>
            <div class="cols compact" style="--n: 2; margin-top: 0.25em">
                <div class="panel">
                    <h3>Core distance κ(x)</h3>
                    <p class="small">
                        The distance to the k-th neighbour: small where stars
                        are packed, large where they are sparse. Two of them
                        give the <strong>mutual reachability</strong> distance
                        max(κ(x), κ(y), d(x, y)) that HDBSCAN* clusters on.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>The k-NN graph</h3>
                    <p class="small">
                        Keep the edges, not just the distance: every point
                        joined to its k neighbours. That graph is the object
                        <strong>UMAP</strong> and <strong>EVoC</strong> embed;
                        clustering becomes a question about a graph.
                    </p>
                </div>
            </div>
            <aside class="notes">
                (~4 min) The pivot slide; everything after this is a caller of
                it, so do not rush. KNN is normally taught as a supervised
                classifier, but here we steal its one useful idea: the
                neighbourhood. Two things fall out of it. First, the
                k-th-neighbour distance κ(x) is a free local density estimate;
                that is DBSCAN's minPts test and HDBSCAN*'s mutual reachability.
                Second, the neighbour lists are a graph, and that graph is what
                UMAP and EVoC embed. Point at the figure: same k, two very
                different radii, and that difference is the density signal.
            </aside>
        </section>

        <!-- 12a · KNN step 1 -->
        <section class="dense">
            <div class="eyebrow">
                Tool two &middot; KNN &middot; step 1 of 4
            </div>
            <h2>Pick k, then measure every distance</h2>
            <div
                class="figure"
                style="width: 74%; max-height: 52vh; margin: 0.3em auto 0"
            >
                <img
                    :src="asset('knn_step1_distances.png')"
                    alt="One point joined to every other point in the set, showing the exhaustive distance computation"
                    style="width: 100%; height: auto; display: block"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 820px; margin-inline: auto"
            >
                The only knob is <strong>k</strong>. Naively every query costs
                O(N); a k-d or ball tree brings it down to about O(log N), which
                is why the primitive is cheap enough to call constantly.
            </p>
            <aside class="notes">
                (~40 s) Start deliberately dumb: to find the nearest neighbours
                you first measure everything. Say out loud that nobody
                implements it this way (the trees are what make it practical)
                but the definition really is this simple, and that simplicity is
                why four later methods can lean on it.
            </aside>
        </section>

        <!-- 12b · KNN step 2 -->
        <section class="dense">
            <div class="eyebrow">
                Tool two &middot; KNN &middot; step 2 of 4
            </div>
            <h2>Keep the k smallest; that is the neighbourhood</h2>
            <div
                class="figure"
                style="width: 74%; max-height: 52vh; margin: 0.3em auto 0"
            >
                <img
                    :src="asset('knn_step2_neighbourhood.png')"
                    alt="The same point now joined only to its six nearest neighbours"
                    style="width: 100%; height: auto; display: block"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 820px; margin-inline: auto"
            >
                One list per star. Everything that follows (density,
                reachability, graphs) is a different <em>use</em> of this one
                list.
            </p>
            <aside class="notes">
                (~30 s) This is the whole algorithm. The interesting part is not
                computing the list, it is what you decide the list means, and
                the next two slides give two completely different answers.
            </aside>
        </section>

        <!-- 12c · KNN step 3 -->
        <section class="dense">
            <div class="eyebrow">
                Tool two &middot; KNN &middot; step 3 of 4
            </div>
            <h2>Reading 1; the k-th distance is a free density estimate</h2>
            <div
                class="figure"
                style="width: 74%; max-height: 52vh; margin: 0.3em auto 0"
            >
                <img
                    :src="asset('knn_step3_core_distance.png')"
                    alt="Two query points with the same k: the dense one has a small core-distance circle, the sparse one a much larger circle"
                    style="width: 100%; height: auto; display: block"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 880px; margin-inline: auto"
            >
                Same k, different radius. Small &kappa;(x) = packed, large
                &kappa;(x) = sparse. <strong>DBSCAN</strong> thresholds this;
                <strong>HDBSCAN*</strong> builds mutual reachability
                max(&kappa;(x), &kappa;(y), d(x,y)) from it.
            </p>
            <aside class="notes">
                (~50 s) The pivot. You never asked for a density estimate, but
                you got one for free: fix the count and let the radius float,
                and the radius <em>is</em> the density. Point at both circles,
                same six neighbours, radius roughly three times bigger in the
                sparse blob. Every density method in this deck is downstream of
                this one picture.
            </aside>
        </section>

        <!-- 12d · KNN step 4 -->
        <section class="dense">
            <div class="eyebrow">
                Tool two &middot; KNN &middot; step 4 of 4
            </div>
            <h2>Reading 2; keep the edges and you have a graph</h2>
            <div
                class="figure"
                style="width: 74%; max-height: 52vh; margin: 0.3em auto 0"
            >
                <img
                    :src="asset('knn_step4_graph.png')"
                    alt="Every point joined to its nearest neighbours, forming a k-nearest-neighbour graph over the whole dataset"
                    style="width: 100%; height: auto; display: block"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 880px; margin-inline: auto"
            >
                Now clustering is a <strong>graph</strong> question, not a
                geometry question. This object is what
                <strong>UMAP</strong> lays out and what
                <strong>EVoC</strong> clusters.
            </p>
            <aside class="notes">
                (~45 s) Second reading of the same list, and the one that
                carries the back half of the talk. Notice what was thrown away:
                absolute distances. The graph keeps only who-is-near-whom, which
                is exactly the robustness UMAP exploits in high dimensions. One
                primitive, two products, density on the previous slide, topology
                on this one.
            </aside>
        </section>

        <!-- 13 · DBSCAN -->
        <section>
            <div class="eyebrow">Tool three · DBSCAN</div>
            <h2>DBSCAN: follow the density</h2>
            <div class="slide-body">
                <div
                    class="fig-split"
                    style="
                        --cols: 1.3fr 1fr;
                        margin-top: 0.3em;
                        align-items: start;
                    "
                >
                    <div>
                        <p class="small">
                            No centres at all: clusters are high-density
                            regions, grown from neighbour counts instead of
                            nearest-centre distance (Ester et al. 1996). Two
                            parameters, three kinds of point.
                        </p>
                        <ol class="contribs tight small">
                            <li>
                                <strong>Pick</strong> the neighbourhood radius
                                <strong>ε</strong> and the minimum count
                                <strong>minPts</strong>.
                            </li>
                            <li>
                                <strong>Core point:</strong> x is
                                <em>core</em> if its ε-ball holds at least
                                minPts points, counting x itself.
                            </li>
                            <li>
                                <strong>Grow:</strong> a cluster is a core point
                                plus everything density-reachable from it, a
                                chain of cores, each within ε of the next.
                            </li>
                            <li>
                                <strong>Border point:</strong> a non-core point
                                inside a core's ε-ball joins that cluster, but
                                never extends it.
                            </li>
                            <li>
                                <strong>Noise:</strong> reachable from no core →
                                outlier. DBSCAN need not cluster everything.
                            </li>
                        </ol>
                    </div>
                    <div>
                        <div class="figure" style="margin: 0">
                            <img
                                :src="asset('dbscan.gif')"
                                alt="DBSCAN growing clusters: core points with radius epsilon, density-reachable points joining, outliers marked with a cross"
                                style="
                                    width: 100%;
                                    height: auto;
                                    display: block;
                                "
                            />
                        </div>
                        <p class="small muted" style="margin-top: 0.5em">
                            The KNN thread: the core test is still a neighbour
                            count, <strong>k fixed, radius free</strong> in KNN;
                            <strong>radius fixed, k free</strong> in DBSCAN.
                        </p>
                    </div>
                </div>
            </div>
            <aside class="notes">
                (~3 min) Define both knobs out loud before anything else: ε is
                how far you look, minPts is how many you need to find. Then walk
                the animation, a core point recruits its ε-ball, the recruits
                that are themselves core keep the chain going, the ones that are
                not sit on the border and stop it. DBSCAN's two big wins over
                K-means: arbitrary shapes (no centroid, no Gaussian), and
                uniquely so far; it labels outliers instead of forcing every
                point into a cluster. That "some points are noise" option is
                exactly right for chemical tagging, where most field stars
                belong to nothing. Close on the KNN thread: same neighbour
                count, only which of k and the radius is held fixed changes. But
                DBSCAN has its own assumption, coming next.
            </aside>
        </section>

        <!-- 14 · Pick ε and minPts -->
        <section>
            <div class="eyebrow">Tool three · DBSCAN · step 1 of 4</div>
            <h2>Pick ε and minPts</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('dbscan_step1_knobs.png')"
                    alt="Two points, each with a circle of radius epsilon drawn around it: one holds five points and passes the minPts test, the other holds three and fails"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                ε is a radius, so draw it: put a ball of that size on a point
                and count what falls inside. Pass minPts and the point is
                <em>core</em>.
            </p>
            <aside class="notes">
                (~30 s) Two knobs, and both are visible here. ε is how far you
                look; literally the circle. minPts is how many you need to find
                inside it, counting the point itself. Show the pass and the near
                miss together: five clears the bar of four, three does not. Same
                neighbour-counting idea as KNN, only now the radius is fixed and
                the count varies, instead of the other way round.
            </aside>
        </section>

        <!-- 15 · Classify every point: core, border, noise -->
        <section>
            <div class="eyebrow">Tool three · DBSCAN · step 2 of 4</div>
            <h2>Classify every point: core, border, noise</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('dbscan_step2_classify.png')"
                    alt="The same points, now labelled: twelve core points shown with their epsilon-balls, three border triangles inside a core's ball, four noise crosses inside nobody's"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Draw a ball on every core and the three roles read straight off:
                inside your own ball, inside someone else's, or inside nobody's.
            </p>
            <aside class="notes">
                (~30 s) One rule, applied to all nineteen points, and it sorts
                them into three piles. Core: the ball around it is crowded
                enough. Border: it fails the count itself, but it sits inside a
                core's ball, so it joins that cluster. Noise: no core's ball
                reaches it. Only the twelve cores can start or extend anything;
                that is the sentence to land before the next slide.
            </aside>
        </section>

        <!-- 16 · Grow clusters from core points -->
        <section>
            <div class="eyebrow">Tool three · DBSCAN · step 3 of 4</div>
            <h2>Link the cores, and the clusters appear</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('dbscan_step3_grow.png')"
                    alt="Double-headed arrows between core points within epsilon of each other form two connected chains; single-headed arrows run out to the border triangles; noise points sit inside empty dashed balls"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Two cores within ε link <strong>both ways</strong>; a core links
                to a border <strong>one way only</strong>. A cluster is one
                connected chain, and nobody chose how many.
            </p>
            <aside class="notes">
                (~30 s) The arrows are the definition. Between two cores the
                link runs both ways, so the chain keeps going; out to a border
                it runs one way only, which is exactly why a border joins a
                cluster but can never grow it. Follow a chain and you have a
                cluster: two of them here, and nothing chose that number; it is
                however many connected chains the data happens to contain. That
                is the answer to K-means asking you for K.
            </aside>
        </section>

        <!-- 17 · Read off shapes and outliers -->
        <section>
            <div class="eyebrow">Tool three · DBSCAN · step 4 of 4</div>
            <h2>Read off shapes and outliers</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('dbscan_step4_result.png')"
                    alt="DBSCAN on 300 points of interleaved half-moons plus uniform noise: both crescents recovered whole, with the scattered points marked as noise crosses"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Same two knobs, hundreds of stars, and these are the half-moons
                K-means cut down the middle on the limits slide.
            </p>
            <aside class="notes">
                (~30 s) Same ε, same minPts, no new machinery, just more points.
                Two wins over K-means, and both are on this one picture.
                Arbitrary shapes: these are the interleaved crescents that
                K-means sliced in half back on the limits slide, and a chain of
                ε-balls follows a curve as happily as a blob. And noise is an
                allowed answer; the scattered points are left out rather than
                forced into whichever cluster is nearest. For chemical tagging
                that second one is the point: most field stars belong to no
                cluster at all.
            </aside>
        </section>

        <!-- 18 · No free lunch -->
        <section>
            <div class="eyebrow">The toolbox, compared</div>
            <h2>No free lunch</h2>
            <p class="small">
                K-means and DBSCAN each win on their own toy data and lose on
                the others, a small change flips the ranking, and averaged over
                all problems no algorithm is universally best. KNN is not a
                clusterer: it is the <strong>primitive</strong> underneath
                DBSCAN and everything after.
            </p>
            <div class="cols" style="--n: 3; margin-top: 0.4em">
                <div class="panel">
                    <h3>K-means</h3>
                    <p class="small">
                        Fast, simple. Assumes spherical, similar-size clusters;
                        needs K; Euclidean → linear boundaries.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>KNN</h3>
                    <p class="small">
                        The primitive. "Who's close?", no clusters, just
                        neighbours; feeds DBSCAN, HDBSCAN*, UMAP, EVoC.
                    </p>
                </div>
                <div class="panel">
                    <h3>DBSCAN</h3>
                    <p class="small">
                        Arbitrary shape + outliers. Assumes one density
                        threshold works for every cluster.
                    </p>
                </div>
            </div>
            <p class="stat center" style="margin-top: 0.4em">
                Know the assumptions
            </p>
            <p class="small center muted">
                The job is to match them to your data, and to check with
                validation.
            </p>
            <aside class="notes">
                (~2 min) The no-free-lunch beat (the BigDataLDN through-line).
                No model is "bad"; each encodes assumptions, and the failure is
                applying a model outside them. KNN is the odd one out by design:
                it makes no clusters at all, it just measures locality, which is
                why it shows up inside every method that follows. This is the
                philosophical core the whole school wants you to internalise,
                and it leads straight into the density section.
            </aside>
        </section>

        <!-- 19 · Validation — internal -->
        <section>
            <div class="eyebrow">Validation · internal scores</div>
            <h2>How do you know the clusters are real?</h2>
            <p class="small">
                <strong>Internal</strong> validation uses only the data and the
                labels the algorithm just produced, no answer key. Three scores,
                three different questions, and a sharp caveat on each.
            </p>
            <div class="cols" style="--n: 3; margin-top: 0.4em">
                <div class="panel">
                    <h3>SSE</h3>
                    <p class="small">
                        Summed squared distance from every object to its own
                        cluster mean. It falls with every extra cluster, so it
                        cannot choose K, but across random restarts the run with
                        the lowest SSE is the one to keep.
                    </p>
                </div>
                <div class="panel">
                    <h3>Silhouette</h3>
                    <p class="small">
                        Mean within-cluster distance against mean
                        between-cluster distance, averaged over every object: −1
                        to 1, where 1 is cleanly separated and 0 is completely
                        overlapping. Its peak over K suggests K; trust it only
                        when the peak is clear.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Dip test</h3>
                    <p class="small">
                        Run this <em>first</em>: is the distribution multipeaked
                        at all? The null hypothesis is a single peak, so a small
                        p-value means real structure. Silhouette only means
                        something once this passes (Hartigan &amp; Hartigan
                        1985).
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.5em">
                Silhouette gets a slide of its own next. Homogeneity and ARI
                need true labels; that is <em>external</em> validation, two
                slides on.
            </p>
            <aside class="notes">
                (~3 min) The chapter's §13.2.2, in the order you should actually
                use them: (1) is the data multipeaked at all, the dip test, and
                if it fails, stop; (2) how many groups, the silhouette, and only
                if it has a clear peak, because on unstructured data it returns
                essentially random values; (3) which run to keep, the lowest SSE
                across restarts, since SSE falls monotonically with K and can
                never pick K for you. Say the honest bottom line the chapter
                states: these scores are aids, not oracles, and the
                interpretation has to be informed by domain knowledge. Then flag
                the trap on the last line, homogeneity (Rosenberg &amp;
                Hirschberg 2007, and the score the chapter itself uses to
                compare K-means, GMM and DBSCAN) and ARI get quoted as if they
                were internal scores, but both compare against truth, which is a
                completely different kind of check. The luxury of this problem
                is that we have one, which is the next slide.
            </aside>
        </section>

        <!-- 19b · Silhouette, in detail -->
        <section class="denser">
            <div class="eyebrow">Validation · internal · silhouette</div>
            <h2>Silhouette: is this point closer to home than to next door?</h2>
            <div class="figure" style="margin: 0.15em auto 0; max-width: 1020px">
                <img
                    :src="asset('silhouette_explained.png')"
                    alt="Left: one point with dashed rings showing its mean distance to its own cluster and to the nearest other cluster. Right: a silhouette plot with every point's score sorted inside each cluster, the overall mean marked, and the Kaufman and Rousseeuw thresholds."
                    style="width: auto; max-width: 100%; max-height: 330px; height: auto; display: block; margin: 0 auto"
                />
            </div>
            <div class="cols compact" style="--n: 3; margin-top: 0.4em">
                <div class="panel">
                    <h3>Per point, not per cluster</h3>
                    <p class="small">
                        <strong>a</strong> is the mean distance to its own
                        cluster, <strong>b</strong> the mean distance to the
                        nearest cluster it is <em>not</em> in. The score is
                        their gap over the larger of the two, so it lands
                        between &minus;1 and 1 whatever the units.
                    </p>
                </div>
                <div class="panel">
                    <h3>Read the plot, not the average</h3>
                    <p class="small">
                        One mean hides everything. The plot shows each cluster's
                        own spread, and a cluster whose bars run short or
                        negative is the one to distrust, even when the headline
                        number looks healthy.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Where it misleads</h3>
                    <p class="small">
                        It rewards <strong>round, separated</strong> blobs, so
                        it marks down exactly the shapes DBSCAN exists to find,
                        and it needs every pairwise distance. In high dimension
                        the distances converge and the whole scale compresses.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.45em">
                Kaufman &amp; Rousseeuw's published reading:
                <strong>0.71+</strong> strong, <strong>0.51</strong> reasonable,
                <strong>0.26</strong> weak and possibly artificial, below that,
                no substantial structure.
            </p>
            <aside class="notes">
                (~2 min) Walk the left panel first: pick one star, measure the
                average distance to its own cluster (a), then to the nearest
                cluster it does not belong to (b). If b is much bigger than a
                the point is comfortably home and the score approaches 1; if
                they are equal it sits on the border at 0; if a is bigger the
                point is closer to the neighbours than to its own label and the
                score goes negative. The figure's numbers are computed, not
                drawn: a = 0.90, b = 3.34, so s = 0.73 for that point.
                Then the right panel, and this is the part people skip: the
                overall mean here is 0.65, but the value of the silhouette is
                the <em>shape</em> of the plot. Three clusters sit at 0.67,
                0.67 and 0.59; had one been at 0.15 the mean would still look
                respectable while one cluster was junk. Land the two caveats.
                First, the thresholds are stricter than people assume, 0.65 is
                only "reasonable", not "strong". Second, it is a convexity
                score: it prefers round separated blobs, so a low silhouette on
                a half-moon or a filament means the metric disagrees with the
                shape, not that the cluster is fake. That is precisely why the
                dip test comes first and why this deck does not let silhouette
                pick the winner on its own.
            </aside>
        </section>

        <!-- 20 · Validation — external truth -->
        <section>
            <div class="eyebrow">Validation · external truth</div>
            <h2>The luxury of this problem: ground truth</h2>
            <p class="small">
                Internal scores can only say a grouping is tight and separated.
                <strong>External</strong> validation asks the harder question;
                is it the <em>right</em> grouping? and that needs labels the
                clustering never saw.
            </p>
            <div class="cols" style="--n: 3; margin-top: 0.4em">
                <div class="panel">
                    <h3>Where it comes from</h3>
                    <p class="small">
                        Chemistry is not the only fingerprint: stars born
                        together also share position, parallax, proper motion
                        and radial velocity. Gaia kinematics label membership
                        without ever touching the 16 abundances we clustered on.
                    </p>
                </div>
                <div class="panel">
                    <h3>What it buys</h3>
                    <p class="small">
                        <strong>Recall</strong>, did we recover the real
                        members? <strong>Precision</strong>, are the ones we
                        claim real? And homogeneity and ARI, the label-invariant
                        scores, finally mean something.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Why it matters</h3>
                    <p class="small">
                        Most unsupervised problems (brain imaging, say) have no
                        answer key, so the argument never ends. Chemical tagging
                        has one, and that is what turns it into a benchmark: the
                        results later in this talk are measured, not hoped for.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.5em">
                Internal asks <em>is this grouping self-consistent?</em>,
                external asks <em>is this grouping true?</em>
            </p>
            <aside class="notes">
                (~2 min) The pivot from "trust the score" to "check against
                reality", and the contrast with the previous slide is the point:
                internal scores never leave the data you clustered, external
                validation brings in an answer from outside it. Stress the
                independence, the labels come from six-dimensional phase space,
                not from the 16 abundances, so scoring chemistry against
                kinematics is not circular. That is why the benchmark at the end
                means something: we do not have to argue about whose silhouette
                is better, we can quote recall and precision. Keep this in your
                back pocket; it returns at the benchmark slides.
            </aside>
        </section>

        <!-- 20b · Homogeneity, in detail -->
        <section class="denser">
            <div class="eyebrow">Validation · external · homogeneity</div>
            <h2>Homogeneity: does every cluster hold just one kind of star?</h2>
            <div class="figure" style="margin: 0.15em auto 0; max-width: 1020px">
                <img
                    :src="asset('homogeneity_explained.png')"
                    alt="Three clusterings of the same labelled points: one recovers the truth and scores 1.00 on both measures; one shatters each class into three and still scores homogeneity 1.00 while completeness falls to 0.50; one merges two classes and scores completeness 1.00 with homogeneity 0.58."
                    style="width: auto; max-width: 100%; max-height: 330px; height: auto; display: block; margin: 0 auto"
                />
            </div>
            <div class="cols compact" style="--n: 3; margin-top: 0.4em">
                <div class="panel">
                    <h3>Purity, one cluster at a time</h3>
                    <p class="small">
                        Look inside a cluster and ask how mixed the true labels
                        are. All one class is zero uncertainty and scores 1.
                        Formally <strong>h = 1 &minus; H(C|K)/H(C)</strong>, the
                        class entropy left once the cluster is known.
                    </p>
                </div>
                <div class="panel">
                    <h3>It cannot be read alone</h3>
                    <p class="small">
                        Give every star its own cluster and homogeneity is a
                        perfect 1.00, having learned nothing. That is the middle
                        panel. <strong>Completeness</strong> asks the mirror
                        question, and V-measure is the harmonic mean of the two.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Why this deck uses it</h3>
                    <p class="small">
                        It needs no cluster-to-class matching and does not care
                        how the labels are numbered, so it compares runs that
                        found different numbers of clusters. But it moves with
                        the number of clusters, so quote it with a chance level.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.45em">
                Cluster count fixed by the algorithm, not by you, so always
                report homogeneity <em>and</em> completeness; one without the
                other is a number you can game.
            </p>
            <aside class="notes">
                (~2 min) The figure is one labelled set clustered three ways,
                and the scores are computed by scikit-learn, not asserted.
                Start left: recovering the truth gives 1.00 and 1.00, the
                uninteresting case. The middle panel is the one that matters:
                every true class cut into three still scores homogeneity 1.00,
                because every cluster is still pure, while completeness
                collapses to 0.50. Purity is free if you are allowed to cut
                finely enough; in the limit, one star per cluster scores a
                perfect 1.00 and has told you nothing. The right panel is the
                opposite failure, two classes lumped together: completeness
                1.00, homogeneity 0.58. So homogeneity alone is not a result,
                it is half of one. Connect it forward twice: this is the pair
                behind the V-measure numbers in the benchmark table, and it is
                why we quote a chance level there, a random partition into
                similarly sized groups already scores well above zero. It is
                also the score the chapter itself uses to compare K-means, GMM
                and DBSCAN, so the audience will meet it again.
            </aside>
        </section>

        <!-- 21 · The scale problem -->
        <section>
            <div class="eyebrow">The scale problem</div>
            <h2>One density threshold can't see everything</h2>
            <p class="small">
                DBSCAN commits to a single density threshold ε. But stellar
                structure lives across a
                <strong>huge range of densities</strong>, and one ε can only
                ever pick one point on that axis.
            </p>
            <p class="small muted" style="margin: 0.45em 0 0.1em">
                Sparse → dense, all in the same survey:
            </p>
            <div class="scale-bar"></div>
            <div class="scale-ticks">
                <span>tidal tails</span>
                <span>moving groups</span>
                <span>the field disc</span>
                <span>open clusters</span>
                <span>globular cores</span>
            </div>
            <div class="cols" style="--n: 3; margin-top: 0.55em">
                <div class="panel">
                    <h3>The symptom</h3>
                    <p class="small">
                        Loosen ε until the tidal tail hangs together, and the
                        disc fuses into one blob. Tighten it until the globular
                        splits cleanly, and the tail is all noise.
                    </p>
                </div>
                <div class="panel">
                    <h3>The half-fix</h3>
                    <p class="small">
                        <strong>OPTICS</strong> (Ankerst et al. 1999) refuses to
                        choose: it plots the density profile at every scale and
                        lets the analyst read the valleys. Honest, but still a
                        picture a human has to interpret, one dataset at a time.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>What we want</h3>
                    <p class="small">
                        A <em>hierarchy over all density scales</em>, with every
                        cluster extracted at the scale where it actually exists
                        and the reading-off done for us.
                    </p>
                </div>
            </div>
            <aside class="notes">
                (~2 min) The transition slide. Say it plainly: DBSCAN's one
                assumption (a single global density threshold) is the first
                thing real data breaks. Walk the bar left to right and name a
                real object at each end so they feel the dynamic range: a
                dissolving tidal tail is orders of magnitude sparser than a
                globular core, and no single ε serves both. Then OPTICS as the
                honest halfway house; it shows you the whole profile instead of
                picking for you, but a human still has to read the picture, and
                that does not scale to a survey. The next two slides automate
                exactly that reading: HDBSCAN* builds the hierarchy, PLSCAN
                selects from it.
            </aside>
        </section>

        <!-- 22 · HDBSCAN* -->
        <section class="dense">
            <div class="eyebrow">Density at all scales · HDBSCAN*</div>
            <h2>HDBSCAN*: a hierarchy over every density</h2>
            <p class="small">
                HDBSCAN* (Campello et al. 2013, 2015; McInnes &amp; Healy 2017)
                runs DBSCAN at
                <strong>every threshold at once</strong> and keeps the whole
                tree of answers.
            </p>
            <div class="cols compact" style="--n: 2; margin-top: 0.25em">
                <div class="panel">
                    <h3>The machinery</h3>
                    <p class="small">
                        <strong>Mutual reachability.</strong> d<sub>mut</sub> =
                        max(κ<sub>i</sub>, κ<sub>j</sub>, d<sub>ij</sub>), with
                        κ the distance to the k-th neighbour. Gaps shorter than
                        κ are levelled up, pushing sparse points away.
                    </p>
                    <p class="small">
                        <strong>Build, then read.</strong> Single-linkage on
                        d<sub>mut</sub> gives one branch per structure over the
                        density range it survives; slicing at one density
                        recovers DBSCAN.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>The knob you still pick</h3>
                    <p class="small">
                        Two survive: <strong>k</strong> (behind κ) and
                        <strong>min_cluster_size m<sub>c</sub></strong
                        >, how many stars a branch needs to count as a cluster
                        rather than a wiggle.
                    </p>
                    <p class="small">
                        m<sub>c</sub> is a <strong>smoothing knob</strong>:
                        raise it and shallow peaks disappear, as if the density
                        were blurred (Bot et al. 2025, Fig. 1, a 2-D toy cloud,
                        not stars). A small m<sub>c</sub> lets a handful of
                        stars count as a cluster; a large one prunes that peak
                        away and the same stars become noise.
                    </p>
                </div>
            </div>
            <div
                class="fig-split"
                style="
                    --cols: 1.54fr 1.29fr;
                    max-width: 640px;
                    margin: 0.25em auto 0;
                "
            >
                <div class="figure" style="width: 100%; box-sizing: border-box">
                    <img
                        :src="asset('mutual_reachability.png')"
                        alt="Two points, each ringed by its core-distance circle, with the straight-line distance between them and the max-of-three definition of mutual reachability"
                        style="width: 100%; height: auto"
                    />
                </div>
                <div class="figure" style="width: 100%; box-sizing: border-box">
                    <img
                        :src="asset('hdbscan_density.gif')"
                        alt="Animation sweeping the density level: clusters appear, grow and merge as the threshold drops"
                        style="width: 100%; height: auto"
                    />
                </div>
            </div>
            <aside class="notes">
                (~3 min) Two ideas, in this order. First mutual reachability,
                pointing at the left figure: taking the max of the two core
                distances and the true distance means any gap smaller than a
                point's own core distance gets levelled up to it. Sparse points
                end up far from everything, which is exactly what stops
                single-linkage chaining through noise. Second, the sweep on the
                right; narrate it as "this is DBSCAN, at every ε,
                simultaneously". The tree records the level at which each clump
                is born and the level at which it merges away, and the branch
                that survives the longest span is the one you want. Then set up
                the next slide: the density threshold is gone, but m_c is not,
                and m_c is a genuine judgement call with no data-driven answer.
                (The right-hand panel is an animation; if you are presenting
                from the PDF export, talk through the sweep instead of pointing
                at it.)
            </aside>
        </section>

        <!-- 23 · Core distance κ(x) -->
        <section class="denser">
            <div class="eyebrow">Tool four · HDBSCAN* · step 1 of 4</div>
            <h2>Core distance κ(x)</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 740px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('knn_core_distance.png')"
                    alt="Two query points, each ringed by its core-distance circle"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                The distance to the k-th neighbour, the local density probe you
                already met in KNN.
            </p>
            <aside class="notes">
                (~30 s) Reuse the KNN slide's picture deliberately: κ is the
                same quantity. Dense regions give small κ, sparse regions large
                κ.
            </aside>
        </section>

        <!-- 24 · Mutual reachability d_mut -->
        <section>
            <div class="eyebrow">Tool four · HDBSCAN* · step 2 of 4</div>
            <h2>Mutual reachability d_mut</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('mutual_reachability.png')"
                    alt="Two points with core-distance circles and the max-of-three definition"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                d_mut = max(κ_i, κ_j, d_ij): gaps shorter than κ are levelled
                up, sparse points are pushed apart.
            </p>
            <aside class="notes">
                (~30 s) The key move. Taking the max flattens dense regions and
                separates sparse points; this is what stops single-linkage from
                chaining through noise.
            </aside>
        </section>

        <!-- 25 · Build one tree over all densities -->
        <section>
            <div class="eyebrow">Tool four · HDBSCAN* · step 3 of 4</div>
            <h2>Build one tree over all densities</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('hdbscan_step3_tree.png')"
                    alt="A single-linkage dendrogram with a horizontal density cut line"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Single-linkage on d_mut gives one branch per structure, slicing
                at any density recovers DBSCAN.
            </p>
            <aside class="notes">
                (~30 s) One tree contains every DBSCAN run. The horizontal cut
                is one density threshold; the tree keeps all of them.
            </aside>
        </section>

        <!-- 26 · Keep the longest-lived branches -->
        <section>
            <div class="eyebrow">Tool four · HDBSCAN* · step 4 of 4</div>
            <h2>Keep the longest-lived branches</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('hdbscan_step4_sweep.png')"
                    alt="Clusters visible at one density level of the sweep"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Run every density threshold at once; the branches that survive
                the longest span are your clusters.
            </p>
            <aside class="notes">
                (~30 s) The sweep is DBSCAN at every ε simultaneously.
                Longest-lived branch wins; that is the only selection rule, and
                m_c is the knob left over for next slide.
            </aside>
        </section>

        <!-- 27 · PLSCAN -->
        <section class="denser">
            <div class="eyebrow">Persistence · PLSCAN</div>
            <h2>PLSCAN: drop the min-cluster-size knob</h2>
            <p class="small">
                <strong>PLSCAN</strong>, Persistent Leaves Spatial Clustering
                for Applications with Noise (Bot, McInnes &amp; Aerts 2025)
                keeps HDBSCAN*'s hierarchy but stops picking a scale. Instead it
                measures
                <strong
                    >how long each cluster survives across all of them</strong
                >.
            </p>
            <div class="cols compact" style="--n: 3; margin-top: 0.45em">
                <div class="panel">
                    <h3>1 · Scale-space</h3>
                    <p class="small">
                        Raising m<sub>c</sub> never moves a merge; it only
                        prunes branches too small to count. So one condensed
                        tree already contains every HDBSCAN* run at every
                        m<sub>c</sub>.
                    </p>
                </div>
                <div class="panel">
                    <h3>2 · Leaf tree</h3>
                    <p class="small">
                        For each leaf cluster, record the interval of m<sub
                            >c</sub
                        >
                        over which it is a leaf. Its length, s<sub>max</sub> −
                        s<sub>min</sub>, is that cluster's
                        <strong>persistence</strong>, one bar of the barcode.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>3 · Persistence trace</h3>
                    <p class="small">
                        Sum the persistence of every leaf alive at each
                        m<sub>c</sub>. Peaks in that trace are the most stable
                        scales, and PLSCAN returns one
                        <strong>layer of detail</strong> per peak.
                    </p>
                </div>
            </div>
            <p class="small muted" style="margin: 0.45em 0 0">
                Formally, 0-dimensional persistent homology on a purpose-built
                metric space (their App. A). Practically: no m<sub>c</sub> to
                guess. <strong>k survives</strong>, which is the next slide's
                business.
            </p>
            <div
                class="fig-split"
                style="
                    --cols: 1.42fr 2.39fr;
                    max-width: 900px;
                    margin: 0.4em auto 0;
                "
            >
                <div class="figure" style="width: 100%; box-sizing: border-box">
                    <img
                        :src="asset('plscan_persistence.gif')"
                        alt="Sweeping min-cluster-size across a three-peak density profile; each peak's bar grows for as long as that cluster stays alive"
                        style="width: 100%; height: auto"
                    />
                </div>
                <div class="figure" style="width: 100%; box-sizing: border-box">
                    <img
                        :src="asset('plscan_barcode.png')"
                        alt="A three-peak density profile beside its persistence barcode: one bar per cluster, its length the range of min-cluster-size over which that cluster exists"
                        style="width: 100%; height: auto"
                    />
                </div>
            </div>
            <aside class="notes">
                (~3 min) The newest material in the talk, so go slowly and lean
                on the pictures. Start from last slide's complaint: m_c is a
                smoothing knob with no right answer. PLSCAN's move is to refuse
                the question, sweep m_c over its whole range and ask which
                clumps are still there at the end. Run the left animation and
                watch peaks drop out as the threshold rises. The right figure is
                the same information redrawn: one bar per cluster, bar length =
                the span of m_c over which it survived = its persistence. Here
                the densest peak persists to 202 and the shallowest only to 70;
                say that a bar this short would be discarded as a fluctuation on
                real data, but note the toy figure only draws the three
                survivors, not the stubs. Be honest about the caveat in the
                muted line: PLSCAN is not parameter-free, k is still there. What
                it removes is the parameter you had no principled way to set.
            </aside>
        </section>

        <!-- 28 · Find the density peaks -->
        <section>
            <div class="eyebrow">Tool five · PLSCAN · step 1 of 3</div>
            <h2>Find the density peaks</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('plscan_step1_density.png')"
                    alt="A three-peak density profile with each peak marked"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Every density peak is a candidate cluster, no min_cluster_size
                to guess.
            </p>
            <aside class="notes">
                (~30 s) The starting point: peaks in the density profile. PLSCAN
                refuses to pick a scale and instead considers all of them.
            </aside>
        </section>

        <!-- 29 · Measure each cluster's persistence -->
        <section>
            <div class="eyebrow">Tool five · PLSCAN · step 2 of 3</div>
            <h2>Measure each cluster's persistence</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('plscan_step2_barcode.png')"
                    alt="One horizontal bar per cluster, length equal to its persistence"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                One bar per cluster: its length is how long that cluster
                survives as min_cluster_size sweeps.
            </p>
            <aside class="notes">
                (~30 s) Persistence = lifetime in min_cluster_size. A bar this
                long means the cluster is still there across a wide range of
                scales.
            </aside>
        </section>

        <!-- 30 · Keep the long bars, discard the stubs -->
        <section>
            <div class="eyebrow">Tool five · PLSCAN · step 3 of 3</div>
            <h2>Keep the long bars, discard the stubs</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('plscan_step3_read.png')"
                    alt="Barcode with long bars coloured and short bars greyed, a cut line between"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Long bars are stable clusters; short bars are fluctuations, the
                data, not you, draws the line.
            </p>
            <aside class="notes">
                (~30 s) The long/short split replaces the m_c knob you used to
                set by hand. Caveat: k still survives; PLSCAN is not
                parameter-free.
            </aside>
        </section>

        <!-- 31 · PLSCAN — results & segue -->
        <section class="dense">
            <div class="eyebrow">Persistence · PLSCAN</div>
            <h2>More stable, less sensitive, and a familiar name</h2>
            <div class="slide-body">
                <p class="small">
                    Bot et al. benchmark PLSCAN against HDBSCAN* (excess-of-mass
                    selection) across a suite of real-world datasets, at the
                    conventional default k = 4. Medians over those datasets:
                </p>
                <div class="kpis" style="--n: 4; margin-top: 0.45em">
                    <div class="panel center">
                        <div class="stat">0.66</div>
                        <div class="stat-label">
                            median ARI<br />HDBSCAN*: 0.57
                        </div>
                    </div>
                    <div class="panel center">
                        <div class="stat">0.85</div>
                        <div class="stat-label">
                            V-measure<br />HDBSCAN*: 0.79
                        </div>
                    </div>
                    <div class="panel center">
                        <div class="stat pink">0.74</div>
                        <div class="stat-label">
                            homogeneity<br />HDBSCAN*: 0.82
                        </div>
                    </div>
                    <div class="panel center">
                        <div class="stat pink">0.76</div>
                        <div class="stat-label">
                            non-noise fraction<br />HDBSCAN*: 0.88
                        </div>
                    </div>
                </div>
                <div class="cols compact" style="--n: 3; margin-top: 0.5em">
                    <div class="panel">
                        <h3>Less sensitive to k</h3>
                        <p class="small">
                            HDBSCAN*'s labels swing about at low k. PLSCAN
                            returns much the same clustering across the tested
                            range.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Affordable</h3>
                        <p class="small">
                            Run-times competitive with K-means++ on
                            low-dimensional data; at high dimension it scales
                            like HDBSCAN*, which is the price of using a space
                            tree at all.
                        </p>
                    </div>
                    <div class="panel flip">
                        <h3>The trade-off</h3>
                        <p class="small">
                            The pink numbers are the bill: PLSCAN's clusters are
                            more complete (0.93 vs 0.89) but less pure, and it
                            sends more stars to the noise bin. A better ARI is
                            not "better everywhere", and every score is computed
                            <strong>only over non-noise points</strong>, so
                            PLSCAN's 0.66 is measured on the 76% it keeps
                            against HDBSCAN*'s 88%.
                        </p>
                    </div>
                </div>
                <p class="small center muted" style="margin-top: 0.5em">
                    Middle author: <strong>Leland McInnes</strong>, also behind
                    UMAP and the HDBSCAN* implementation everyone actually runs.
                    Which is where we go next.
                </p>
            </div>
            <aside class="notes">
                (~2 min) Read the four numbers as one sentence, not four: PLSCAN
                wins on the summary scores because it recovers more of each true
                cluster, and it pays for that with lower purity and more stars
                thrown to noise. That is the deck's through-line again, no free
                lunch, only different assumptions. Do not oversell it: this is
                one benchmark suite at one default k, on general-purpose
                datasets rather than on abundances, and the authors say
                themselves that whether density maxima are the "true" clusters
                depends on the use case. Then the segue: the middle author is
                the person behind UMAP and behind the HDBSCAN* library, which is
                why EVoC later fuses all three; one research programme, not
                three coincidences. And now we change strategy entirely: stop
                clustering in the native space, and reshape the space first.
            </aside>
        </section>

        <!-- 32 · t-SNE — project to 2-D -->
        <section class="dense">
            <div class="eyebrow">Embeddings · t-SNE</div>
            <h2>t-SNE: project to 2-D, then look</h2>
            <p class="small">
                Everything so far clustered in the native abundance space. Kos
                et al. (2017) change the question:
                <strong>reshape the space first</strong>. We are excellent at
                spotting a clump in ≤ 3-D and hopeless in 13-D, so t-distributed
                stochastic neighbour embedding compresses the chemistry into a
                picture, and clustering becomes drawing a line round what you
                can already see.
            </p>
            <div
                class="fig-split"
                style="--cols: 1fr 1fr; align-items: start; margin-top: 0.5em"
            >
                <div class="panel">
                    <h3>Their recipe, on GALAH</h3>
                    <p class="small">
                        <strong>13 abundances, Manhattan distance.</strong>
                        Summing |Δ| instead of squaring it, so one bad line does
                        less damage than under Euclidean.
                    </p>
                    <p class="small">
                        <strong>Weights ∝ 1 / cluster scatter.</strong> Ba and
                        K, whose scatter is mostly noise, drop to 0.25; the
                        clean Fe, Ti, Cr and Cu rise to 2.0, a stand-in for
                        per-star errors they never measured.
                    </p>
                    <p class="small">
                        <strong>Perplexity.</strong> How many neighbours count
                        as "local", next slide.
                    </p>
                    <p class="small">
                        <strong
                            >One sky region per cluster, 30–45° in
                            radius.</strong
                        >
                        9,408 stars within 40° of the Pleiades, say, never the
                        whole survey at once.
                    </p>
                </div>
                <div class="figure" style="width: 100%; box-sizing: border-box">
                    <img
                        :src="asset('tsne.gif')"
                        alt="A t-SNE embedding converging from a scattered layout into four separated clusters"
                        style="width: 100%; height: auto"
                    />
                </div>
            </div>
            <aside class="notes">
                (~3 min) Flag the change of strategy first; it is the hinge of
                the lecture. Up to now we asked an algorithm to find clusters in
                high dimension; from here we compress to two and hand the
                problem to the audience's visual cortex. Then the three
                practical choices that are easy to skim past and matter
                enormously. Manhattan distance for outlier robustness. The
                weights, which are the inverse of the observed cluster scatter,
                Ba and K scatter mostly because they are hard to measure, so
                they get down-weighted, while Fe, Ti, Cr and Cu are trusted; it
                is a poor man's error bar, and Kos et al. say the map falls
                apart without it. And above all the region cut: they never ran
                one all-sky t-SNE, every map is a 30–45° cone around one known
                cluster. Promise them the region cut comes back to bite us in
                the benchmark.
            </aside>
        </section>

        <!-- 33 · t-SNE — how it works -->
        <section>
            <div class="eyebrow">Embeddings · t-SNE</div>
            <h2>t-SNE, under the hood</h2>
            <p class="small">
                t-SNE never tries to preserve distances. It turns both spaces
                into probability distributions over "who is whose neighbour",
                then matches them.
            </p>
            <div
                class="cols compact"
                style="
                    --n: 3;
                    grid-template-columns: 1.3fr 1fr 1fr;
                    margin-top: 0.5em;
                "
            >
                <div class="panel">
                    <h3>The mechanics</h3>
                    <ul class="dotlist small">
                        <li>
                            High-D: Gaussian p<sub>ij</sub>, one bandwidth σ<sub
                                >i</sub
                            >
                            per star
                        </li>
                        <li>
                            Low-D: heavy-tailed Student-t q<sub>ij</sub> on the
                            map
                        </li>
                        <li>
                            Minimise <strong>KL(P‖Q)</strong> by gradient
                            descent
                        </li>
                        <li>Barnes-Hut brings the cost to O(N log N)</li>
                        <li>Random start, so every run differs</li>
                        <li>
                            Standard practice: keep the lowest-KL map of many
                            though Kos et al. report their own repeated runs
                            differed only by a random rotation
                        </li>
                    </ul>
                </div>
                <div class="panel">
                    <h3>Perplexity: the knob</h3>
                    <p class="small">
                        Each σ<sub>i</sub> is tuned until star i's neighbour
                        distribution has the <strong>perplexity</strong> you
                        asked for, an effective neighbour count, typically 5–50.
                        It plays the role k plays in KNN.
                    </p>
                    <p class="small">
                        Low → fine local clumps. High → coarse global shape.
                        Same data, different picture, so always quote it with
                        the map.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>The crowding problem</h3>
                    <p class="small">
                        A 13-D ball has far more room at middling distances than
                        any 2-D disc. Match Gaussians in both spaces and all
                        those middling neighbours pile inwards into one blob.
                    </p>
                    <p class="small">
                        The Student-t's <strong>heavy tail</strong> is the fix:
                        a decent q<sub>ij</sub> is reachable at a long map
                        distance, so points can spread and gaps open between
                        clumps.
                    </p>
                </div>
            </div>
            <aside class="notes">
                (~2 min) The one mechanics slide; keep it to two takeaways.
                Takeaway one is perplexity: it is the local/global dial, roughly
                "how many neighbours count as near", and a map is meaningless
                without it quoted. Give the 5–50 range and point people at the
                Wattenberg demo that Kos et al. cite. Takeaway two is the
                crowding problem, because it explains the "t" in the name: if
                the map used a Gaussian too, there simply is not enough room in
                two dimensions for everything that sits at moderate distance in
                thirteen, and the whole plot collapses inward. The heavy tail
                buys that room back, and the visible gaps between clumps are its
                doing. Mention the random initialisation in passing; it is why
                two runs look different, and why anyone showing you a t-SNE map
                owes you the perplexity. Save the harder criticisms (cluster
                sizes and inter-cluster distances mean nothing) for the limits
                slide.
            </aside>
        </section>

        <!-- 34 · Model who is whose neighbour, in high-D -->
        <section>
            <div class="eyebrow">Embeddings · t-SNE · step 1 of 3</div>
            <h2>Model who is whose neighbour, in high-D</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('tsne_step1_highd.png')"
                    alt="Two Gaussian similarity curves, narrow and wide"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Each star gets a Gaussian over its neighbours, tuned to a fixed
                perplexity, an effective neighbour count.
            </p>
            <aside class="notes">
                (~30 s) High-D similarities p_ij. Perplexity plays the role k
                played in KNN: how many neighbours count as near.
            </aside>
        </section>

        <!-- 35 · Model the same in 2-D, with a heavy tail -->
        <section>
            <div class="eyebrow">Embeddings · t-SNE · step 2 of 3</div>
            <h2>Model the same in 2-D, with a heavy tail</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('tsne_step2_lowd.png')"
                    alt="A Gaussian curve beside a heavy-tailed Student-t curve"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                A Student-t in the map, the heavy tail lets points spread
                instead of crowding into one blob.
            </p>
            <aside class="notes">
                (~30 s) The crowding problem and its fix. The 't' in t-SNE is
                this Student-t; the visible gaps between clumps are its doing.
            </aside>
        </section>

        <!-- 36 · Match the two, by gradient descent -->
        <section>
            <div class="eyebrow">Embeddings · t-SNE · step 3 of 3</div>
            <h2>Match the two, by gradient descent</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('tsne_step3_descent.png')"
                    alt="Final t-SNE embedding with four separated clusters"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Minimise KL(P‖Q): similar stars are pulled together, dissimilar
                ones pushed apart.
            </p>
            <aside class="notes">
                (~30 s) Gradient descent matches the two distributions. Random
                start, so every run differs, quote the perplexity with the map.
            </aside>
        </section>

        <!-- 37 · t-SNE results -->
        <section>
            <div class="eyebrow">t-SNE results · Kos et al. 2017</div>
            <h2>It recovered clusters, and found new Pleiades</h2>
            <p class="small muted" style="margin: 0.1em 0 0.3em">
                GALAH (pilot + main survey) + K2-HERMES,
                <strong>13 abundances</strong>, 9408 stars in a 40° radius
                around the Pleiades. t-SNE draws the map; the groups are then
                drawn on it <em>by eye</em>.
            </p>
            <div
                class="fig-split"
                style="--cols: 1fr 1.4fr; align-items: start"
            >
                <div>
                    <div
                        class="figure"
                        style="aspect-ratio: 1 / 1; max-height: 50vh"
                    >
                        <img
                            :src="asset('kos_pleiades_tsne.png')"
                            alt="Kos et al. 2017 Fig. 2, members panel: t-SNE map of 9408 stars, field grey, Pleiades members red in two tight groups A and B, the two new members green and numbered 1 and 2 inside a blue polygon"
                            style="
                                width: 100%;
                                height: 100%;
                                object-fit: contain;
                            "
                        />
                    </div>
                    <p
                        class="small muted"
                        style="margin: 0.2em 0 0; font-size: 0.55em"
                    >
                        Kos et al. 2017, Fig. 2 (members panel). Red: known
                        members; green 1, 2: the new ones; blue: the hand-drawn
                        group.
                    </p>
                </div>
                <div>
                    <ul class="checklist small">
                        <li>
                            <strong>7 of 9</strong> clusters recovered from
                            chemistry alone
                        </li>
                        <li>Nine observed: six globular, three open</li>
                        <li>
                            <strong>Two new Pleiades members</strong>,
                            kinematically confirmed; one had been a supercluster
                            candidate before; the other is new to the literature
                        </li>
                        <li>One of them 6° out, a full tidal radius</li>
                        <li>
                            Sub-groups of a cluster land together on the map
                        </li>
                    </ul>
                    <ul class="crosslist small">
                        <li>
                            <strong>47 Tuc</strong>, no chemical group at all
                            (the paper can only guess why)
                        </li>
                        <li>
                            <strong>NGC 2516</strong>, only an edge observed, so
                            no verdict
                        </li>
                        <li>Some field-star contamination in the groups</li>
                        <li>
                            <em
                                >“a large fraction of the stars are
                                untaggable”</em
                            >
                        </li>
                    </ul>
                </div>
            </div>
            <aside class="notes">
                (~3 min) This is the proof the approach works on real stars, say
                the numbers out loud, then spend equal time on the right-hand
                column. Two of the nine failed for completely different reasons:
                47 Tuc simply has no chemical group in the map, while NGC 2516
                was barely observed, so it is a data failure, not a method
                failure. Worth flagging that the Pleiades sub-groups A and B the
                paper sees in The Cannon abundances could not be confirmed with
                SME; the hierarchy claim is about the projection preserving
                structure, not about that particular split being real. Close on
                the paper's own verdict: with 13 elements a large fraction of
                stars are untaggable. Chemical tagging is hard which is exactly
                why we benchmark it. If asked what "recovered" means: a majority
                of members land in one group, not all of them, 17 of the 27
                known Pleiades members sit in groups A and B (the Fig. 2
                caption). The new-member logic is a two-step funnel: chemistry
                cuts 9408 stars to about 30 candidates, then radial velocity,
                proper motion and distance cut those 30 to 2. The paper's chance
                estimate for a coincidence is 2 × 30 / 9400 ≈ 0.006 stars. One
                of the two (star 2) was already a supercluster candidate that
                their own membership cut had missed; star 1 has no previous link
                to the Pleiades.
            </aside>
        </section>

        <!-- 38 · t-SNE limits -->
        <section>
            <div class="eyebrow">t-SNE limits</div>
            <h2>t-SNE embeds; it does not cluster</h2>
            <div class="cols" style="--n: 2; margin-top: 0.3em">
                <div class="panel">
                    <h3>What it won't tell you</h3>
                    <ul class="crosslist small">
                        <li>
                            <strong>Distance</strong>, gaps between clumps are
                            arbitrary
                        </li>
                        <li>
                            <strong>Density</strong>, area on the map isn't
                            volume
                        </li>
                        <li>
                            <strong>Size</strong>, sparse groups spread, dense
                            ones shrink
                        </li>
                        <li>
                            <strong>Stability</strong>, rerun it and the clumps
                            move
                        </li>
                        <li>
                            <strong>Speed</strong>, O(N²), even with Barnes-Hut
                        </li>
                    </ul>
                </div>
                <div class="panel flip">
                    <h3>Where do labels come from?</h3>
                    <p class="small">
                        Not from t-SNE. It hands you a <em>picture</em>, Kos et
                        al. drew the polygons by hand. To automate that you bolt
                        a clusterer onto the map:
                        <strong>HDBSCAN* on the 2-D embedding</strong>.
                    </p>
                    <ul class="crosslist small">
                        <li>Two algorithms, two sets of knobs</li>
                        <li>It clusters the picture, not the data</li>
                        <li>So it inherits every distortion on the left</li>
                    </ul>
                </div>
            </div>
            <p class="small center" style="margin-top: 0.4em">
                Read a t-SNE map for <strong>neighbourhoods</strong>, never
                distance, size or density.
            </p>
            <p class="small center muted" style="margin-top: 0.15em">
                Can one method embed <em>and</em> cluster, and pick the scale
                itself?
            </p>
            <aside class="notes">
                (~2 min) This is the slide to slow down on: the single most
                common mistake with t-SNE is reading it as a map with a scale.
                Say it plainly; t-SNE embeds, it does not cluster, and neither
                the distances nor the densities in the picture are trustworthy.
                The mitigation is the two-step pipeline, which works but
                inherits every weakness of both halves, because HDBSCAN* is
                clustering the distortion, not the abundances. The last line
                sets up the two requirements EVoC will eventually meet: one
                method, and a scale it chooses itself.
            </aside>
        </section>

        <!-- 39 · UMAP -->
        <section class="dense">
            <div class="eyebrow">Embeddings · UMAP</div>
            <h2>UMAP: a graph embedding with global structure</h2>
            <p class="small muted" style="margin: 0.1em 0 0.3em">
                Uniform Manifold Approximation and Projection (McInnes, Healy
                &amp; Melville 2018) skips the N² pairs entirely: it builds a
                <strong>weighted neighbour graph</strong>, then draws that
                graph. Hold on to the graph; EVoC reuses it, and so does its
                author's other work.
            </p>
            <div
                class="fig-split"
                style="--cols: 1.6fr 1fr; align-items: start"
            >
                <div>
                    <ol class="contribs small tight">
                        <li>
                            <strong>kNN graph</strong>, join each star to its
                            <strong>n_neighbors</strong> neighbours.
                        </li>
                        <li>
                            <strong>Fuzzy edges</strong>, rescale each point's
                            distances by its own nearest-neighbour distance, so
                            a weight is the <em>probability an edge exists</em>,
                            not a distance.
                        </li>
                        <li>
                            <strong>Layout</strong>, attract along edges, repel
                            sampled non-neighbours, minimising cross-entropy
                            against that graph (t-SNE minimises KL).
                        </li>
                    </ol>
                    <ul class="dotlist small" style="margin-top: 0.35em">
                        <li>
                            <strong>n_neighbors</strong>, how local the graph is
                            (cf. perplexity)
                        </li>
                        <li>
                            <strong>min_dist</strong>, how tightly the layout
                            may pack points
                        </li>
                    </ul>
                </div>
                <div>
                    <div class="figure">
                        <img
                            :src="asset('umap.gif')"
                            alt="UMAP embedding converging from a spectral layout to four separated clusters"
                            style="width: 100%; height: auto"
                        />
                    </div>
                    <ul class="checklist small" style="margin-top: 0.4em">
                        <li>Near-linear in N</li>
                        <li>Keeps global layout</li>
                        <li>Still a projector</li>
                        <li>HDBSCAN* still after</li>
                    </ul>
                </div>
            </div>
            <aside class="notes">
                (~3 min) Walk the three numbered steps slowly; this is the one
                place in the talk where a graph, not a distance matrix, is the
                object being optimised, and everything after it depends on that.
                The fuzzy step is the one to dwell on: because each point's
                distances are rescaled by its own nearest neighbour, a sparse
                region and a dense region get comparable edge weights, which is
                how UMAP avoids t-SNE's density distortion. Then the payoff: it
                is near-linear rather than quadratic, and unlike t-SNE the
                layout of clusters relative to each other means something. But
                it is still only a projector; HDBSCAN* still has to run
                afterwards. Land the last line hard: EVoC starts from this same
                graph.
            </aside>
        </section>

        <!-- 40 · Build the k-nearest-neighbour graph -->
        <section>
            <div class="eyebrow">Embeddings · UMAP · step 1 of 3</div>
            <h2>Build the k-nearest-neighbour graph</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('umap_step1_graph.png')"
                    alt="Points joined by edges to their nearest neighbours"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Join each star to its k neighbours; the graph is the object
                everything after this is drawn on.
            </p>
            <aside class="notes">
                (~30 s) The same graph EVoC will reuse. Unlike t-SNE, UMAP never
                builds the N² distance matrix; it only keeps these edges.
            </aside>
        </section>

        <!-- 41 · Make the edges fuzzy -->
        <section>
            <div class="eyebrow">Embeddings · UMAP · step 2 of 3</div>
            <h2>Make the edges fuzzy</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('umap_step2_fuzzy.png')"
                    alt="Edges with thickness proportional to their fuzzy weight"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Rescale each edge by local density, so a weight means
                'probability this edge exists'.
            </p>
            <aside class="notes">
                (~30 s) The fuzzy step is the one to dwell on: rescaling by each
                point's own nearest neighbour makes sparse and dense regions
                comparable, how UMAP avoids t-SNE's density distortion.
            </aside>
        </section>

        <!-- 42 · Lay the graph out -->
        <section>
            <div class="eyebrow">Embeddings · UMAP · step 3 of 3</div>
            <h2>Lay the graph out</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 780px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('umap_step3_layout.png')"
                    alt="Final UMAP embedding with four separated clusters"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 760px;
                    margin-inline: auto;
                "
            >
                Attract along edges, repel sampled non-neighbours, the clusters
                separate.
            </p>
            <aside class="notes">
                (~30 s) The layout keeps global structure, unlike t-SNE. But it
                is still a projector; HDBSCAN* runs afterwards. Bridge to EVoC.
            </aside>
        </section>

        <!-- 43 · EVoC — the fusion -->
        <section>
            <div class="eyebrow">The fusion · EVoC</div>
            <h2>EVoC: embed and cluster in one pass</h2>
            <p class="small">
                <strong>EVōC</strong>, Embedding Vector Oriented Clustering
                (Tutte Institute; said “evoke”) is not a projector. It is the
                three ideas we have just built, fused into a single fit that
                hands back <strong>labels, not coordinates</strong>.
            </p>
            <div class="cols" style="--n: 3; margin-top: 0.45em">
                <div class="panel">
                    <h3>from UMAP</h3>
                    <p class="small">
                        The same kNN graph on the raw 16 abundances, laid out by
                        the same node-embedding machinery, but into 4–15
                        dimensions of its own, never a 2-D picture.
                    </p>
                </div>
                <div class="panel">
                    <h3>from HDBSCAN*</h3>
                    <p class="small">
                        Mutual-reachability MST and condensed tree on that
                        embedding: density clustering, with a noise label for
                        the field stars that belong to nothing.
                    </p>
                </div>
                <div class="panel">
                    <h3>from PLSCAN</h3>
                    <p class="small">
                        Persistence over min_cluster_size scores every layer of
                        the hierarchy and returns the most persistent one. No
                        scale to guess.
                    </p>
                </div>
            </div>
            <div class="cols" style="--n: 1; margin-top: 0.4em">
                <div class="panel flip">
                    <p class="small" style="margin: 0">
                        <strong>The caveat to hold on to:</strong> EVoC works in
                        <strong>cosine</strong>
                        geometry. Our benchmark L2-normalises every abundance
                        vector so Euclidean ≈ cosine; that is what keeps the
                        three methods comparable, and it turns out to be the
                        single biggest precision lever in the whole study.
                    </p>
                </div>
            </div>
            <aside class="notes">
                (~4 min) This is the payoff of the whole tour, so make the three
                columns explicit: nothing in EVoC is new; column one is the UMAP
                graph, column two is HDBSCAN*, column three is PLSCAN. What is
                new is that they share one fit, so the embedding is built for
                the clustering rather than for your eyes. Correct the obvious
                misreading before it happens: EVoC does not cluster the raw 16-D
                vectors, it clusters its own internal node embedding, the point
                is that there is no lossy 2-D detour, not that the geometry is
                untouched. Then the caveat: cosine is its native metric, and
                rather than call that unfair we normalise the vectors so all
                three methods see the same geometry. Next slide opens the hood.
            </aside>
        </section>

        <!-- 44 · EVoC — the pipeline -->
        <section>
            <div class="eyebrow">The fusion · EVoC</div>
            <h2>EVoC, under the hood</h2>
            <ol class="contribs small tight" style="margin-top: 0.2em">
                <li>
                    <strong>kNN graph</strong>, n_neighbors neighbours per star,
                    on all 16 abundances, in cosine geometry.
                </li>
                <li>
                    <strong>Node embedding</strong>, UMAP-style layout of that
                    graph, label-propagation init, into 4–15-D.
                </li>
                <li>
                    <strong>Mutual-reachability MST</strong>, Borůvka on the
                    embedding: HDBSCAN*'s own step, unchanged.
                </li>
                <li>
                    <strong>Condensed tree</strong>, the density hierarchy at
                    every scale, exactly as HDBSCAN* builds it.
                </li>
                <li>
                    <strong>Persistence</strong>, a min_cluster_size barcode
                    scores each layer; the most persistent one wins.
                </li>
            </ol>
            <p class="small muted" style="margin-top: 0.3em">
                Out:
                <strong
                    >labels + membership strengths + persistence scores + the
                    layer hierarchy</strong
                >
                not one flat partition.
            </p>
            <div class="center" style="margin-top: 0.3em">
                <div
                    class="figure"
                    style="aspect-ratio: 1444 / 443; width: 90%"
                >
                    <img
                        :src="asset('evoc_pipeline.png')"
                        alt="EVoC pipeline: kNN graph, node embedding, cluster labels, persistence per layer"
                        style="width: 100%; height: 100%; object-fit: contain"
                    />
                </div>
            </div>
            <aside class="notes">
                (~3 min) Read the five steps as a map back onto the tour: 1–2 is
                UMAP's graph machinery, 3–4 is HDBSCAN*, 5 is PLSCAN. Two
                details are worth pausing on. First, the density clustering does
                not run on the 16-D abundances; it runs on the node embedding,
                which EVoC sizes itself (four dimensions at our n_neighbors=15).
                Second, step 5 is literally the PLSCAN barcode: it sweeps
                min_cluster_size, scores each resulting layer by total
                persistence, and the winning layer becomes the labels you get
                back. Point at the four panels as you go; the last one is the
                persistence score per layer.
            </aside>
        </section>

        <!-- 44a · EVoC step 1 -->
        <section class="dense">
            <div class="eyebrow">
                The fusion &middot; EVoC &middot; step 1 of 4
            </div>
            <h2>Build the kNN graph on the raw abundances</h2>
            <div
                class="figure"
                style="width: 72%; max-height: 52vh; margin: 0.3em auto 0"
            >
                <img
                    :src="asset('evoc_step1_graph.png')"
                    alt="A k-nearest-neighbour graph drawn over the full dataset in its original feature space"
                    style="width: 100%; height: auto; display: block"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 880px; margin-inline: auto"
            >
                Borrowed wholesale from <strong>UMAP</strong>, in cosine
                geometry, on all 16 dimensions. Note what does
                <em>not</em> happen here: no 2-D picture is ever made.
            </p>
            <aside class="notes">
                (~40 s) Step 1 is the KNN slide again, unchanged, say that
                explicitly, it is reassuring. The one thing to stress is the
                metric: cosine, because the signal is the abundance pattern
                rather than its amplitude.
            </aside>
        </section>

        <!-- 44b · EVoC step 2 -->
        <section class="dense">
            <div class="eyebrow">
                The fusion &middot; EVoC &middot; step 2 of 4
            </div>
            <h2>Embed the graph, into 4–15-D, not into a picture</h2>
            <div
                class="figure"
                style="width: 72%; max-height: 52vh; margin: 0.3em auto 0"
            >
                <img
                    :src="asset('evoc_step2_embed.png')"
                    alt="The graph laid out as a node embedding, with the groups now separated"
                    style="width: 100%; height: auto; display: block"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 880px; margin-inline: auto"
            >
                Sized for the <em>clusterer</em>, not for your eyes, four
                dimensions at our n_neighbors&nbsp;=&nbsp;15. Drawn in 2-D here
                only so it can be put on a slide.
            </p>
            <aside class="notes">
                (~45 s) The key distinction in the whole method, and the one the
                audience will get wrong if you let them: this is not a
                visualisation. t-SNE and UMAP squash to 2-D because a human has
                to look; EVoC embeds to whatever dimension clusters best. Be
                honest that the slide cheats by showing 2-D; the real thing has
                no picture at all.
            </aside>
        </section>

        <!-- 44c · EVoC step 3 -->
        <section class="dense">
            <div class="eyebrow">
                The fusion &middot; EVoC &middot; step 3 of 4
            </div>
            <h2>Run HDBSCAN* on that embedding</h2>
            <div
                class="figure"
                style="width: 72%; max-height: 52vh; margin: 0.3em auto 0"
            >
                <img
                    :src="asset('evoc_step3_cluster.png')"
                    alt="The embedded points coloured by cluster, with sparse points labelled as noise"
                    style="width: 100%; height: auto; display: block"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 880px; margin-inline: auto"
            >
                Mutual-reachability MST, then the condensed tree,
                <strong>HDBSCAN*</strong>'s own steps, unchanged. Field stars
                are allowed to be noise, which is what a cluster search needs.
            </p>
            <aside class="notes">
                (~40 s) Nothing new here at all, and that is the point, steps
                1–3 are two methods they have already seen, bolted together.
                Point at the grey noise points: K-means could never produce
                them, and in a catalogue that is ~97% field stars that ability
                is the whole game.
            </aside>
        </section>

        <!-- 44d · EVoC step 4 -->
        <section class="dense">
            <div class="eyebrow">
                The fusion &middot; EVoC &middot; step 4 of 4
            </div>
            <h2>The new idea, score every layer, keep the best</h2>
            <div
                class="figure"
                style="width: 76%; max-height: 50vh; margin: 0.3em auto 0"
            >
                <img
                    :src="asset('evoc_step4_persistence.png')"
                    alt="Bar chart of mean cluster stability against min_cluster_size, with the winning layer highlighted"
                    style="width: 100%; height: auto; display: block"
                />
            </div>
            <p
                class="small muted center"
                style="margin-top: 0.3em; max-width: 900px; margin-inline: auto"
            >
                Sweep min_cluster_size, score each layer by
                <strong>persistence</strong>, keep the winner. Numbers above the
                bars are clusters found: the smallest setting shatters the data
                into 45 specks and scores badly; the winning layer recovers the
                true structure. <strong>No scale left to guess.</strong>
            </p>
            <aside class="notes">
                (~60 s) This is the only genuinely new component, and it is
                PLSCAN's idea applied to whole layers rather than to individual
                clusters. Walk the bars left to right: tiny min_cluster_size
                fragments everything, and fragmentation is penalised because we
                score the <em>mean</em>
                persistence; a shattered layer is full of short-lived junk. The
                winner is the layer whose typical cluster survives the longest
                span of density thresholds. That is how the method eliminates
                the knob every earlier algorithm made you guess.
            </aside>
        </section>

        <!-- 45 · EVoC — why C-space -->
        <section class="dense">
            <div class="eyebrow">The fusion · EVoC</div>
            <h2>Why this fits chemical space</h2>
            <div class="slide-body">
                <div
                    class="fig-split"
                    style="
                        --cols: 1.3fr 1fr;
                        margin-top: 0.25em;
                        align-items: start;
                    "
                >
                    <div>
                        <ul class="checklist small">
                            <li>
                                <strong>Clusters are rare and tiny.</strong>
                                Tens of members hiding in ~10⁵ field stars; you
                                need a method allowed to call most of the data
                                noise. K-means must give every star a home.
                            </li>
                            <li>
                                <strong
                                    >Groups are elongated, and their densities
                                    differ wildly.</strong
                                >
                                Compact metal-poor globulars and diffuse open
                                clusters in one catalogue: no single ε or
                                min_cluster_size serves both, so persistence
                                picks the layer instead.
                            </li>
                            <li>
                                <strong
                                    >The signal is the abundance pattern, not
                                    its amplitude.</strong
                                >
                                Cosine geometry (EVoC's native metric) is the
                                physically right choice, and it is why we
                                L2-normalise.
                            </li>
                            <li>
                                <strong>No 2-D detour.</strong> The graph is
                                built on all 16 abundances and the clustering
                                runs in EVoC's own 4–15-D embedding, never in a
                                picture.
                            </li>
                        </ul>
                        <p class="small muted" style="margin-top: 0.3em">
                            Its knobs (noise_level, base_min_cluster_size,
                            n_neighbors, n_epochs) are the ones you already met
                            in UMAP and HDBSCAN*.
                        </p>
                    </div>
                    <div class="figure" style="aspect-ratio: 589 / 492">
                        <img
                            :src="asset('cspace_corner.png')"
                            alt="Two of the sixteen abundance dimensions: field stars grey, cluster members magenta in one tight clump"
                            style="
                                width: 100%;
                                height: 100%;
                                object-fit: contain;
                            "
                        />
                    </div>
                </div>
            </div>
            <aside class="notes">
                (~2 min) Bring the corner plot back from the data slide; it is
                the shape of the problem in two of the sixteen dimensions, and
                every bullet on the left is a property you can see in it: a
                handful of red points, elongated, sitting inside a grey cloud
                that is thousands of times larger. Each line pairs one property
                of C-space with the ingredient of EVoC that answers it, so this
                is the slide where the whole narrative arc closes. Then break
                the spell deliberately: this is an argument from assumptions,
                and assumptions are exactly what the talk has been warning
                about. The benchmark is where the story meets reality.
            </aside>
        </section>

        <!-- 45b · Transition — part 1 → part 2. The day's division is the
         speaker's own: 90 min introduction, 90 min problem, 120 min hands-on.
         It lands on the seam the previous slide already set up ("That is the
         benchmark."), so this divider only names the change of mode instead of
         re-arguing part 1. The download check is repeated here because it is
         the first natural pause since the setup slide. -->
        <section class="dense title-slide center">
            <div class="eyebrow">End of part 1 · 90 minutes done</div>
            <h1>The toolbox is full. Now the problem.</h1>
            <p class="subtitle">
                Eight methods, each one fixing the last one's failure, and one
                question still open:
                <strong
                    >does any of it actually recover clusters we already
                    know?</strong
                >
            </p>
            <div class="slide-body">
                <p class="small">
                    Part 2 changes two things. The stars are
                    <strong>real</strong>, with ground truth, so the methods get
                    graded rather than admired. And the judge is
                    <strong>kinematics</strong>, not chemistry.
                </p>
                <p class="small">
                    First, three short slides on how we work together, licence,
                    contributing guide, issues, and the standards that keep a
                    project alive. Then the problem.
                </p>
                <p class="small muted" style="margin-top: 0.4em">
                    Download check · nothing to do, it keeps going. If it has
                    already finished, <code>./run.sh run --fast</code> (~2 min,
                    no GPU) tells you it worked.
                </p>
            </div>
            <aside class="notes">
                (~1 min) This is the hinge of the day: everything so far was the
                toolbox, everything after it is judged. Say it in one breath and
                do not re-argue the previous 90 minutes, the last slide already
                ended on the question ("that is the benchmark") so the only job
                here is to name the change of mode. Two practical things do fit
                naturally at this boundary: ask for hands on the download state
                (it is the first pause since the setup slide) and, if you take a
                break, say exactly when you restart. Then walk into the
                benchmark.
            </aside>
        </section>

        <!-- 45c · Collaboration — the framing slide for the three-slide
         sequence the user asked for at the end of part 1. His own sentence
         ("software is a social contract written in machine language") is the
         heading; the point is that a public repository states its contract in
         the same few places every time, and that today's exercise is a real
         instance of it — the hands-on ends with a pull request to the school
         repo, which CI tests and a maintainer merges. -->
        <section>
            <div class="eyebrow">
                Working together · three slides before part 2
            </div>
            <h2>Software is a social contract written in machine language</h2>
            <p class="small">
                The compiler reads the code. The <strong>contract</strong> is
                what the next person reads: who may use this work, how changes
                are proposed, what is already known. Every public repository
                states it in the same places, and reading them takes minutes.
            </p>
            <div
                class="cols"
                style="--n: 2; margin-top: 0.5em; align-items: start"
            >
                <div class="panel">
                    <h3>What you sign when you contribute</h3>
                    <p class="small">
                        You are joining a project's rules, not just its code:
                        the licence says what may be done with the work, the
                        contributing guide says how they want it done, the
                        issues say what is already in flight. Skipping them is
                        how a good change gets closed unread.
                    </p>
                </div>
                <div class="panel">
                    <h3>What the contract buys you</h3>
                    <p class="small">
                        A change somebody can actually merge: a branch with one
                        reason, a description that says why, tests that pass
                        without your machine, and your name on the work
                        afterwards, the difference between code in a folder and
                        code in a project.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.5em">
                Today's hands-on ends with exactly this: a pull request to the
                school repository, reviewed, tested by CI, merged.
            </p>
            <aside class="notes">
                (~2 min) Introduce it with the sentence on the slide, then make
                it concrete: code is read by machines once and by people
                forever, so what matters is the part that tells other people the
                rules. Two minutes of reading prevents the most common beginner
                mistake in open source, a change nobody asked for, submitted in
                a way nobody asked for. Say plainly that this is not abstract
                for them: at the end of the day their work goes back to the
                school repository as a pull request, the same mechanism as any
                project they will ever contribute to, with the same gates, CI,
                review, merge. Keep it to the two panels and move on; the next
                slide is the practical one.
            </aside>
        </section>

        <!-- 45d · Collaboration — the four pages to read, now one screenshot per
         slide (the user's call: the 4-up row made them too small). Each capture
         is 1490 x 939 and rendered in the theme's native figure mode: natural
         aspect, width auto, capped at 440px tall — 698 x 440 canvas px, about
         2.4x the 4-up row (292 x 178). Two traps to remember: an em cap does
         not bite here (em is ~30px, so a 46em frame went full width at
         1080 x 680 and overflowed by 248px), and the astro theme already caps
         `.dense .figure img` at 60vh with width:auto — and 60vh tracks the real
         browser window, not the 720px canvas, so a `width:100%` image inside a
         frame is letterboxed by object-fit and renders smaller than the frame
         suggests. Explicit px on max-height sidesteps both. The first three are
         astropy's repository and the last is the school repo the students will
         fork, so the sequence reads "a mature project, then yours". The eyebrow
         carries the running order ("n of 4") because these four slides are one
         lesson. -->
        <section class="dense">
            <div class="eyebrow">
                Working together · read the basics · 1 of 4, may I?
            </div>
            <h2>The licence: what you may do with the work</h2>
            <p class="small">
                In this order: <strong>may I</strong> → <strong>how</strong> →
                <strong>is it already known</strong>, then fork. First stop, the
                licence: BSD-3 here, use it, change it, ship it; keep the notice
                and the copyright.
            </p>
            <div style="text-align: center; margin-top: 0.5em">
                <div class="figure" style="padding: 0.3em">
                    <img
                        :src="asset('license_screen.png')"
                        alt="Astropy's LICENSE on GitHub: the BSD-3-Clause summary and full text"
                        style="display: block; max-height: 440px; width: auto"
                    />
                </div>
            </div>
            <p class="caption center">
                astropy/astropy, LICENSE (BSD-3-Clause): permissions,
                limitations, conditions
            </p>
            <aside class="notes">
                (~45 s) The first of the four things to read, and the only one
                with legal weight: can I use this, and what must I keep? Point
                at the three columns (permissions, limitations, conditions) and
                give the BSD-3 shorthand out loud: use it, change it, ship it,
                keep the notice. Then say what no licence would have meant: no
                permission to use it at all, however public the code looks. That
                sentence is worth more than the screenshot.
            </aside>
        </section>

        <!-- 45d2 · Collaboration — the contributing guide, second of four. -->
        <section class="dense">
            <div class="eyebrow">
                Working together · read the basics · 2 of 4, how?
            </div>
            <h2>The contributing guide: how they want it done</h2>
            <p class="small">
                A project that wants help writes down how: what warrants an
                issue, how a change should be shaped, what a pull request has to
                contain, and, in astropy's case, an answer to the doubt.
            </p>
            <div style="text-align: center; margin-top: 0.5em">
                <div class="figure" style="padding: 0.3em">
                    <img
                        :src="asset('contribution_guide.png')"
                        alt="Astropy's CONTRIBUTING.md: reporting issues, contributing code and documentation, and an anti-imposter-syndrome section"
                        style="display: block; max-height: 440px; width: auto"
                    />
                </div>
            </div>
            <p class="caption center">
                astropy/astropy, CONTRIBUTING.md (“We want your help. No,
                really.”)
            </p>
            <aside class="notes">
                (~45 s) This is where a project says how it wants to be helped,
                and it is the page that decides whether your change is welcome
                or wasted: what to report, how to structure a patch, what the
                maintainers will ask for in review. Scroll them to the third
                heading on the screenshot, "anti imposter syndrome reassurance",
                because it answers the objection half the room is silently
                holding: you do not have to be ready, you have to be useful;
                documentation and small fixes count.
            </aside>
        </section>

        <!-- 45d3 · Collaboration — the issues tab, third of four. -->
        <section class="dense">
            <div class="eyebrow">
                Working together · read the basics · 3 of 4, is it already
                known?
            </div>
            <h2>The issues tab: what is already in flight</h2>
            <p class="small">
                Search before you fix. Is it reported? Is somebody already on
                it? Is it deliberate, not a bug? GitHub puts the contributing
                guide on this very page, above the list, for a reason.
            </p>
            <div style="text-align: center; margin-top: 0.5em">
                <div class="figure" style="padding: 0.3em">
                    <img
                        :src="asset('issues_screen.png')"
                        alt="An open issue list, with a banner prompting contributors to read the contributing guidelines before opening an issue"
                        style="display: block; max-height: 440px; width: auto"
                    />
                </div>
            </div>
            <p class="caption center">astropy/astropy, open issues</p>
            <aside class="notes">
                (~45 s) The cheapest ten minutes in open source: the issue list
                tells you what is already known. Read the labels too, "good
                first issue" is how projects mark work they would like a
                newcomer to take, and a bug that is already reported is not a
                bug you should spend the afternoon on. Note the banner GitHub
                shows here, pushing contributors at the contributing guide; that
                is the loop this whole sequence is about.
            </aside>
        </section>

        <!-- 45d4 · Collaboration — the fork, fourth of four. The only one of the
         four screenshots that is not astropy: it is the school repository the
         students will fork today, which is why the notes end by pointing the
         lesson back at their own exercise. -->
        <section class="dense">
            <div class="eyebrow">
                Working together · read the basics · 4 of 4, then fork
            </div>
            <h2>Then fork: your own copy to push to</h2>
            <p class="small">
                Public, Fork, and the licence named in the sidebar. The fork is
                what gives you a branch of your own (and somewhere to push it)
                without asking anyone for permission first.
            </p>
            <div style="text-align: center; margin-top: 0.5em">
                <div class="figure" style="padding: 0.3em">
                    <img
                        :src="asset('fork.png')"
                        alt="The school repository's front page: Public badge, Fork button, and the licence named in the sidebar"
                        style="display: block; max-height: 440px; width: auto"
                    />
                </div>
            </div>
            <p class="caption center">
                github.com/iaa-so-training/iaa-advanced-neural-networks-2026
            </p>
            <aside class="notes">
                (~45 s) Close the sequence on their own repo: this is the page
                they will fork in an hour. Three things to point at, the Public
                badge (you can read and fork it without asking), the Fork button
                (that is where your copy comes from), and the licence in the
                sidebar, which is the same BSD-3 they just saw on astropy. Then
                the sentence that ties the four slides together: every one of
                these pages existed before you arrived, and reading them is the
                whole skill; the code is the easy part.
            </aside>
        </section>

        <!-- 45e · Collaboration — reproducibility and standards, the two things
         that decide whether a contribution survives. Every fact is the workshop
         repository's own: uv with a committed uv.lock and `uv sync --frozen` in
         CI, the Docker path in day_4_clustering/run.sh, pytest (the package has
         229 tests and CI runs the ones needing no real data — the workflow's
         own comment), pyrefly in strict preset and the workflow step that
         labels it "informational — strict mode is not clean yet", and module
         docstrings (cli.py opens with one). details live in docs/
         reproducibility.md and docs/docker.md, which exist in the repo. The
         slide deliberately does not claim a linter the repo does not run: ruff
         is named as what most Python projects reach for, not as this project's
         gate. -->
        <section class="dense">
            <div class="eyebrow">Working together · how it survives</div>
            <h2>Reproducible by default, readable on purpose</h2>
            <div
                class="cols"
                style="--n: 2; margin-top: 0.5em; align-items: start"
            >
                <div class="panel">
                    <h3>Docker + uv · the environment is part of the code</h3>
                    <p class="small">
                        Docker pins the system the code runs in;
                        <code>uv.lock</code> pins every Python package to a
                        version. This workshop does both;
                        <code>./run.sh</code> runs in the container, CI installs
                        with <code>uv sync --frozen</code>, so “it works on my
                        machine” stops being an excuse and becomes a guarantee.
                    </p>
                </div>
                <div class="panel">
                    <h3>Standards · what a reviewer will ask for</h3>
                    <p class="small">
                        Docstrings on public functions, the repo's own CLI opens
                        with one, and type hints, checked by
                        <code>pyrefly</code> in strict preset (still
                        informational, which is honest). Formatting is a tool's
                        job, never a review argument. Plus tests that need no
                        real data: the package's 229 run in CI on every change.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.5em">
                The details live in the repo you already cloned:
                <code>docs/reproducibility.md</code> ·
                <code>docs/docker.md</code>, and this is what your own pull
                request is measured against.
            </p>
            <aside class="notes">
                (~3 min) Two ideas, one line each: the environment is part of
                the code, and style is not a matter of taste. On the first,
                point at what they already did; the container they are running
                and the lockfile behind it are why their run matches the
                person's next to them; `uv sync --frozen` is the same promise in
                CI. On the second, the reviewer's checklist: docstrings so the
                next person knows what a function promises, type hints so the
                mistake surfaces before the test run, formatting done by a tool
                so the diff only contains the change. Be honest about the
                workshop's own state: the type check runs in strict mode but is
                labelled informational because it is not clean yet that is what
                standards look like in real projects, a direction rather than a
                perfect score. Close by pointing at the two docs in their clone;
                they will be graded against the same standard.
            </aside>
        </section>

        <!-- 46 · The benchmark -->
        <section class="denser">
            <div class="eyebrow">The benchmark · this repo</div>
            <h2>Three methods, one kinematic truth</h2>
            <div class="slide-body">
                <p class="small">
                    Reproduce Kos et al. 2017, cluster by cluster, region by
                    region; then extend it:
                    <strong>t-SNE vs UMAP vs EVoC</strong> on APOGEE DR19 + Gaia
                    DR3, over <strong>25 clusters</strong> (18 open, including
                    the Pleiades, and 7 globular; Garcia-Dias et al. 2019 + Kos
                    et al. 2017).
                </p>
                <div class="cols" style="--n: 3; margin-top: 0.5em">
                    <div class="panel">
                        <h3>Pipeline</h3>
                        <p class="small">
                            astraAllStarASPCAP (DR19) FITS → quality cuts (SNR ≥
                            100, ASPCAPFLAG = STARFLAG = 0) → kinematic
                            membership labels → 16-D abundance matrix → embed →
                            cluster → score.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Methods</h3>
                        <p class="small">
                            t-SNE → HDBSCAN* · UMAP → HDBSCAN* · EVoC. The first
                            two only <em>embed</em>, so clustering is a second,
                            separate step; EVoC clusters the 16-D vectors
                            directly and hands back labels.
                        </p>
                    </div>
                    <div class="panel flip">
                        <h3>Score vs the Simbad referee</h3>
                        <p class="small">
                            <strong>recall</strong>, members recovered by the
                            best-matching predicted cluster.
                            <strong>precision</strong>, purity of that cluster.
                            <strong>kNN purity</strong>, share of a member's 10
                            nearest neighbours in the embedding from the same
                            cluster. Ground truth is the
                            <strong>Simbad catalogue</strong> (external
                            literature membership), kinematic labels as
                            fallback.
                        </p>
                    </div>
                </div>
                <p class="small muted center" style="margin-top: 0.55em">
                    C-space: 15 [X/Fe] + [Fe/H], median-centred and unit-scaled,
                    then row-normalised so Euclidean ≡ cosine. Stars need ≥ 8 of
                    16 finite abundances, strict complete-case silently deletes
                    the metal-poor globulars.
                </p>
            </div>
            <aside class="notes">
                (~3 min) This is the repo's whole point: the workshop module.
                t-SNE and UMAP only embed, so HDBSCAN* clusters their 2-D maps;
                EVoC clusters the 16-D vectors internally. Everything is scored
                against the same kinematic ground truth, recall (did we find the
                members?) and precision (were we right?), plus kNN purity as a
                parameter-free cohesion score. Say out loud that kNN purity is
                the automated version of the paper's "draw a polygon round the
                group" test, no hyperparameters, so no way to tune yourself a
                good answer. The footnote matters too: row-normalisation is what
                makes Euclidean distance equal cosine, EVoC's native metric, and
                the imputation rule is what keeps M 15 and M 92 in the sample at
                all.
            </aside>
        </section>

        <!-- 46b · Baseline — re-create the paper -->
        <section class="densest">
            <div class="eyebrow">
                Baseline · re-creating Garcia-Dias et al. 2019
            </div>
            <h2>First: can we separate the clusters from each other?</h2>
            <p class="small">
                The 2019 paper's question, re-run on our DR19 data: take
                <strong>only the known cluster members</strong> (982 stars
                across 25 clusters, one row per star, no field), cluster them,
                and ask how well each star is assigned back to its own cluster,
                scored with the paper's own
                <strong>homogeneity / v-measure / accuracy</strong>.
            </p>
            <div class="cols" style="--n: 2; margin-top: 0.45em">
                <div class="panel">
                    <h3>Abundances only</h3>
                    <table style="font-size: 0.5em; margin-top: 0.3em">
                        <thead>
                            <tr>
                                <th>method</th>
                                <th>homog.</th>
                                <th>v-meas.</th>
                                <th>acc.</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>t-SNE</td>
                                <td>0.56</td>
                                <td>0.58</td>
                                <td>0.44</td>
                            </tr>
                            <tr>
                                <td>UMAP</td>
                                <td><strong>0.58</strong></td>
                                <td>0.58</td>
                                <td>0.39</td>
                            </tr>
                            <tr>
                                <td>EVoC</td>
                                <td>0.42</td>
                                <td>0.47</td>
                                <td>0.40</td>
                            </tr>
                        </tbody>
                    </table>
                    <p class="small muted" style="margin-top: 0.3em">
                        Paper reported 0.85, but on a
                        <em>supervised</em> LDA projection, with membership
                        2&sigma;-clipped in the same abundances (circular). Our
                        honest kinematic membership: 0.42–0.58, mean of 7 seeds.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Abundances + kinematics</h3>
                    <table style="font-size: 0.5em; margin-top: 0.3em">
                        <thead>
                            <tr>
                                <th>method</th>
                                <th>homog.</th>
                                <th>v-meas.</th>
                                <th>acc.</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>t-SNE</td>
                                <td><strong>0.75</strong></td>
                                <td>0.76</td>
                                <td>0.71</td>
                            </tr>
                            <tr>
                                <td>UMAP</td>
                                <td>0.73</td>
                                <td>0.69</td>
                                <td>0.49</td>
                            </tr>
                            <tr>
                                <td>EVoC</td>
                                <td>0.70</td>
                                <td>0.65</td>
                                <td>0.46</td>
                            </tr>
                        </tbody>
                    </table>
                    <p class="small muted" style="margin-top: 0.3em">
                        Add parallax + proper motion + radial velocity:
                        homogeneity jumps to 0.70–0.75. The signal chemistry
                        can't find is in the motion.
                    </p>
                </div>
            </div>
            <p class="stat center" style="margin-top: 0.45em">
                Abundances suggest.<br />Kinematics decide.<br />Giants give
                distance.<br />The main sequence gives age.
            </p>
            <aside class="notes">
                (~3 min) The baseline the whole talk builds on. Re-run the 2019
                question on our data: only cluster members, no field, recover
                which star belongs to which cluster. Two honest corrections to
                the paper's 0.85. First, their best result used LDA, a
                *supervised* projection that already knows the clusters; our
                numbers are fully unsupervised. Second, their membership was
                2σ-clipped in the same abundances it then clustered, circular.
                Our kinematic (Gaia) membership has no such leak, and the honest
                abundances-only number is 0.42–0.58 over seven seeds. Point at
                t-SNE's instability, not just its score: in five of nine row
                orders HDBSCAN* finds just two groups, completeness 0.96,
                homogeneity 0.22, the similar-age open clusters merged into one
                blob, which is the paper's own "indistinguishable pairs". Then
                the flip panel: add parallax, proper motion and radial velocity
                and homogeneity jumps to 0.70–0.75. Chemistry narrows;
                kinematics decide. That is the sentence the rest of the talk
                hangs on, and it lands us on the next slide, pulling clusters
                out of the field, not just apart from each other.
            </aside>
        </section>

        <!-- 47 · Results -->
        <section class="dense">
            <div class="eyebrow">Follow-up · retrieval from the field</div>
            <h2>Chemical tagging is hard, that's the finding</h2>
            <p class="small">
                <strong class="accent"
                    >Recall without precision is a blob, not a
                    discovery.</strong
                >
                In the 30&deg; region sweep UMAP recovers the most members
                (0.56) at precision 0.08; ~92% of its "cluster" is field. t-SNE
                gives recall away and buys the only usable precision: 0.83 on M
                67, which is 10 of that cluster's 269 members. Macro kNN purity:
                0.12 t-SNE, 0.07 UMAP.
            </p>
            <table style="font-size: 0.52em; margin-top: 0.3em">
                <thead>
                    <tr>
                        <th>mode</th>
                        <th>stars</th>
                        <th>t-SNE r / p</th>
                        <th>UMAP r / p</th>
                        <th>EVoC r / p</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>fast, all-sky (abundances)</td>
                        <td>25k</td>
                        <td>0.21 / 0.22</td>
                        <td>0.22 / 0.13</td>
                        <td>0.48 / 0.002</td>
                    </tr>
                    <tr>
                        <td>region sweep (30°), macro over 25 clusters</td>
                        <td>18–63k each</td>
                        <td>0.15 / 0.28</td>
                        <td>0.56 / 0.08</td>
                        <td>0.40 / 0.01</td>
                    </tr>
                    <tr>
                        <td>region, M 67 (30°)</td>
                        <td>40k</td>
                        <td>0.04 / 0.83</td>
                        <td>0.02 / 0.26</td>
                        <td>0.08 / 0.02</td>
                    </tr>
                </tbody>
            </table>
            <div
                class="figure"
                style="
                    display: block;
                    width: calc(100% - 200px);
                    margin: 0.4em 0 0 200px;
                    padding: 0.35em;
                "
            >
                <img
                    :src="asset('benchmark_grid.png')"
                    alt="Three panels: t-SNE and UMAP embeddings with true cluster members coloured and field stars grey, and EVoC's labels drawn on the UMAP canvas"
                    style="display: block; width: 100%; height: auto"
                />
                <div class="caption">
                    Embeddings coloured by true membership, field stars grey.
                    Right panel: EVoC's labels on the UMAP canvas, EVoC clusters
                    in 16-D, it never embeds.
                </div>
            </div>
            <aside class="notes">
                (~3 min) Read the table honestly, row by row. Fast all-sky with
                abundances: everything is mediocre. Region mode (the paper's own
                method) is what moves precision, and t-SNE's 0.83 on M 67 is the
                one number on this slide you could publish; say what it costs:
                that group of 12 stars holds 10 true members, so recall is 0.04,
                or 10 of M 67's 269 members. UMAP's 0.56 recall looks like a win
                until you read across: precision 0.08, so eleven of every twelve
                stars in that "cluster" are field. That is the blob. Then the
                picture: same data three ways, members coloured, field grey.
                Point out that the third panel is not an EVoC projection; EVoC
                has no 2-D output, so its labels are painted onto the UMAP
                canvas. Land it as: chemistry narrows the candidate list,
                kinematics decide.
            </aside>
        </section>

        <!-- 47b · Isochrone + literature -->
        <section>
            <div class="eyebrow">
                Follow-up · isochrone fitting vs the literature
            </div>
            <h2>The same members give back age and distance</h2>
            <p class="small">
                Membership is the input, not the answer. Fit a
                <strong>PARSEC isochrone</strong> to the kinematic members
                (ASteCA + emcee) and read off <strong>age</strong> and
                <strong>distance modulus</strong>, then check against the
                literature (Dias for open clusters, Harris for globulars).
            </p>
            <div class="cols" style="--n: 2; margin-top: 0.45em">
                <div class="panel">
                    <h3>Distance modulus, six seeds, not one</h3>
                    <table style="font-size: 0.52em; margin-top: 0.3em">
                        <thead>
                            <tr>
                                <th>cluster</th>
                                <th>fit dm (6 seeds)</th>
                                <th>lit dm</th>
                                <th>posterior sd</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>M 67 (APOGEE 2-colour)</td>
                                <td>9.51 – 9.81</td>
                                <td>9.54</td>
                                <td>1.1 mag</td>
                            </tr>
                            <tr>
                                <td>M 15 (APOGEE 2-colour)</td>
                                <td>14.35 – 15.14</td>
                                <td>15.04</td>
                                <td>1.1 mag</td>
                            </tr>
                        </tbody>
                    </table>
                    <p class="small muted" style="margin-top: 0.3em">
                        Every best fit lands within 0.7 mag of the literature,
                        but it moves by up to 0.8 mag between seeds, and the
                        posterior is as wide as the flat ±2 mag prior (sd 1.15).
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Age, not constrained yet</h3>
                    <table style="font-size: 0.52em; margin-top: 0.3em">
                        <thead>
                            <tr>
                                <th>cluster</th>
                                <th>fit age, Gyr (6 seeds)</th>
                                <th>lit age</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>M 5 (Gaia CMD)</td>
                                <td>0.9 – 9.8</td>
                                <td>10.5</td>
                            </tr>
                            <tr>
                                <td>M 15 (Gaia CMD)</td>
                                <td>0.03 – 2.2</td>
                                <td>12.5</td>
                            </tr>
                        </tbody>
                    </table>
                    <p class="small muted" style="margin-top: 0.3em">
                        The log-age posterior is as wide as its flat prior (sd ≈
                        1.0 dex): the CMD fit does not constrain age. M 15 also
                        sits on the grid's lowest metallicity node, [M/H] −2.19.
                    </p>
                </div>
            </div>
            <p class="stat center" style="margin-top: 0.45em">
                Distances: in the right place, loosely.<br />Ages: not yet.
            </p>
            <aside class="notes">
                (~2 min) The bridge from "found the members" to "did the
                astrophysics". Kinematic members feed an isochrone fit: ASteCA's
                PARSEC grids + an emcee sampler over (metallicity, log age,
                distance modulus, extinction). Be honest about what it gives.
                Re-run with six seeds, the best-fit distance modulus stays in
                the right place (M 67 9.5–9.8 vs 9.54; M 15 14.4–15.1 vs 15.04).
                The best fit does carry signal: for M 15 it lands a magnitude
                away from the parallax guess of 16.0 that centres the prior, and
                toward the literature. But the chain learns almost nothing
                beyond that point: the posterior is as wide as the flat ±2 mag
                prior (sd 1.1 against 1.15), and over the whole M 67 chain the
                log-probability only moves between −0.45 and −0.11. Treat the
                distance as a sanity check, not a measurement. Ages are not
                constrained at all: for M 67, M 5 and M 15 the log-age posterior
                is as wide as its flat prior (sd ≈ 1.0 dex over 4 Myr, 12.6
                Gyr), so the best fit wanders: M 5 from 0.9 to 9.8 Gyr between
                seeds, M 15 from 0.03 to 2.2 Gyr against 12.5. M 15 also sits on
                the grid's lowest metallicity node ([M/H] −2.19; our literature
                table has [Fe/H] −2.22), so it needs a wider grid as well. The
                next step is a likelihood that actually discriminates (the
                current one is a normalised distance between 0 and 1) before any
                age is quoted. If time, point them at the notebook §7–§10 for
                the HR diagrams and the full literature table.
            </aside>
        </section>

        <!-- 47c-a · Why not just PCA? -->
        <section>
            <div class="eyebrow">
                Dimensionality reduction &middot; the obvious baseline
            </div>
            <h2>Why not just run PCA on the spectrum?</h2>
            <p class="small">
                The spectrum is <strong>8575 pixels</strong> of
                highly-correlated flux. The obvious move is
                <strong>PCA</strong>, project onto the top eigenvectors of the
                pixel covariance. Linear, fast, no labels. Let&rsquo;s see how
                far it gets.
            </p>
            <div class="cols" style="--n: 2; margin-top: 0.4em">
                <div class="panel">
                    <h3>What PCA captures</h3>
                    <ul class="small">
                        <li>
                            the <strong>continuum</strong> (temperature,
                            reddening), the first few PCs
                        </li>
                        <li>global line-to-continuum balance</li>
                        <li>
                            the <strong>largest variance</strong>, not the most
                            informative directions
                        </li>
                    </ul>
                </div>
                <div class="panel flip">
                    <h3>What PCA misses</h3>
                    <ul class="small">
                        <li>
                            each element&rsquo;s
                            <strong>narrow line windows</strong> (low pixel
                            variance)
                        </li>
                        <li>
                            <strong>blends</strong>, overlapping lines from
                            several species
                        </li>
                        <li>
                            the
                            <strong>nonlinear</strong> chemistry&rarr;spectrum
                            mapping
                        </li>
                    </ul>
                </div>
            </div>
            <p class="small muted center" style="margin-top: 0.4em">
                Variance is not chemistry: the continuum dominates the pixel
                variance but carries almost no element-ratio information.
            </p>
            <aside class="notes">
                (~1.5 min) Frame the honest baseline. PCA is the textbook
                answer: 8575 correlated pixels, project to a few eigenvectors.
                It is linear and maximises *variance*. But chemical tagging
                needs *element-ratio* information, and the element ratios live
                in narrow, low-variance absorption lines. PCA spends its budget
                on the continuum (temperature, reddening) because that is where
                the variance is. So PCA is the right baseline and the wrong tool
                the question is whether a learned, nonlinear compression can
                spend its budget on chemistry instead of the continuum.
            </aside>
        </section>

        <!-- 47c-b · The model itself -->
        <section>
            <div class="eyebrow">
                Self-supervised &middot; masked autoencoder
            </div>
            <h2>The masked spectral autoencoder</h2>
            <p class="small" style="margin-bottom: 0.15em; font-size: 0.6em">
                Everything so far consumed <strong>16 abundances</strong>. This
                model consumes the <strong>8575-pixel spectrum</strong> itself
                and learns its own 256 numbers, with
                <strong>no labels at all</strong>.
            </p>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 990px;
                    margin: 0.1em auto 0;
                "
            >
                <img
                    :src="asset('mae_arch.png')"
                    alt="Architecture: 8575-pixel spectrum, five stride-2 conv blocks narrowing to 64 channels, global pool and linear layer to a 256-d latent, then a mirrored transposed-conv decoder back to 8575 pixels"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.25em;
                    font-size: 0.52em;
                    max-width: 900px;
                    margin-inline: auto;
                "
            >
                <strong>8.6M parameters.</strong> Encoder: five
                <code>Conv1d</code> blocks, stride 2, each halving the length
                and widening the channels; then global-average pool and one
                linear layer to <strong>z &isin; &#8477;<sup>256</sup></strong
                >. Decoder mirrors it.
                <strong>z is the only thing we keep.</strong>
            </p>
            <aside class="notes">
                (~3 min) The architecture slide the section was missing, draw it
                out loud, left to right, because everything after this refers
                back to it. The input is one star's spectrum: 8575 flux values,
                standardised to zero mean and unit sigma per star so brightness
                cannot be a feature. Then five convolutional blocks, each stride
                2, so the sequence halves (8575, 4288, 2144, 1072, 536, 268)
                while the channels go the other way: 1024, 512, 256, 128, 64.
                That is the usual convnet trade: lose resolution, gain
                abstraction. Then the step people miss: global-average pool over
                what's left, so the latent does not depend on position, and one
                linear layer down to 256 numbers. The decoder is the mirror
                image, transposed convolutions back up to 8575. Say the
                punchline plainly: at inference we throw the decoder away. The
                decoder exists only to create the training pressure; z is the
                product. And 8.6M parameters is small; this trains in under an
                hour on one GPU, which matters for the hands-on.
            </aside>
        </section>

        <!-- 47c-b1 · Step 1 — mask -->
        <section>
            <div class="eyebrow">
                Self-supervised &middot; masked AE &middot; step 1 of 4
            </div>
            <h2>Hide contiguous blocks of the spectrum</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 900px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('mae_step1_mask.png')"
                    alt="A real DR19 spectrum with two shaded 200-pixel windows hidden from the model"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 820px;
                    margin-inline: auto;
                "
            >
                A real DR19 spectrum. ~50% of it is hidden in
                <strong>200-pixel blocks</strong>, the model is handed the blue
                curve with the shaded windows zeroed out.
            </p>
            <aside class="notes">
                (~45 s) This is a real star, not a cartoon; one of the 39,945
                spectra we actually trained on. Two knobs and both matter: mask
                ratio 50%, block size 200 pixels. Point at a shaded window: the
                model gets zeros there. It has to produce the missing curve from
                everything else it can see. Next slide says why the blocks have
                to be contiguous.
            </aside>
        </section>

        <!-- 47c-b2 · Why blocks, not pixels -->
        <section>
            <div class="eyebrow">
                Self-supervised &middot; masked AE &middot; the design choice
            </div>
            <h2>Why blocks, and not random pixels?</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 980px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('mae_why_blocks.png')"
                    alt="Two panels: with random pixels hidden, linear interpolation reconstructs the spectrum almost perfectly; with one contiguous block hidden, interpolation flatlines across it"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 860px;
                    margin-inline: auto;
                "
            >
                Hide <em>random pixels</em> and plain linear interpolation
                already wins (MSE&nbsp;0.009), the task teaches nothing. Hide
                <em>one block</em> and interpolation flatlines across it
                (MSE&nbsp;0.057, <strong>6&times; worse</strong>).
            </p>
            <aside class="notes">
                (~1.5 min) A slide worth the time, because it is the one design
                decision that makes or breaks a masked model and it generalises
                far beyond spectra. On the left I hid half the pixels at random
                and reconstructed them with nothing but straight-line
                interpolation between the survivors, look how well it does; the
                red dots sit on the blue curve. A model trained on that
                objective learns "average your neighbours", which is a statement
                about sampling, not about stars. On the right I hid one
                contiguous block: interpolation has nothing local left to lean
                on and draws a straight line through real absorption features,
                six times worse. To fill THAT you need to know which lines
                belong there and how deep they are, the iron, the temperature,
                the blends. Same principle as BERT masking whole words and MAE
                masking image patches: the hole must be bigger than the
                correlation length, or the task is trivial.
            </aside>
        </section>

        <!-- 47c-b3 · Step 2 — encode -->
        <section>
            <div class="eyebrow">
                Self-supervised &middot; masked AE &middot; step 2 of 4
            </div>
            <h2>Squeeze 8575 pixels into 256 numbers</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 900px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('mae_step2_encode.png')"
                    alt="Bar chart on a log scale: sequence length falling 8575, 4288, 2144, 1072, 536, 268 while channel count rises 1024 to 64"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 840px;
                    margin-inline: auto;
                "
            >
                Each stride-2 convolution halves the length. After five, a
                <strong>64&times;268</strong> map is pooled over wavelength and
                projected to <strong>z (256-d)</strong>, a
                <strong>33&times;</strong> compression.
            </p>
            <aside class="notes">
                (~1 min) Read the bars as the compression story: 8575 down to
                268 positions, while each position gets richer, 1024 channels at
                the top, 64 at the bottom. The convolutions are local, so early
                layers see individual line profiles and later layers see
                relationships between regions of the spectrum. Then the
                global-average pool: collapse the 268 positions entirely, so the
                latent describes the star rather than a location on the
                detector. 8575 numbers in, 256 out, 33 times smaller, and the
                compression is exactly what forces the model to decide what
                matters.
            </aside>
        </section>

        <!-- 47c-b4 · Step 3 — reconstruct -->
        <section>
            <div class="eyebrow">
                Self-supervised &middot; masked AE &middot; step 3 of 4
            </div>
            <h2>Predict the pixels it was never shown</h2>
            <div
                class="figure"
                style="
                    display: block;
                    width: 100%;
                    max-width: 900px;
                    margin: 0 auto;
                "
            >
                <img
                    :src="asset('mae_step3_recon.png')"
                    alt="The trained model's actual reconstruction, in green, drawn inside the two hidden windows and tracking the true absorption lines"
                    style="width: 100%; height: auto"
                />
            </div>
            <p
                class="small muted center"
                style="
                    margin-top: 0.35em;
                    max-width: 860px;
                    margin-inline: auto;
                "
            >
                The green curve is
                <strong>our trained model's actual output</strong> inside the
                hidden windows. It puts the absorption lines back, in the right
                places, at roughly the right depths.
                <strong>The loss is the MSE there and nowhere else.</strong>
            </p>
            <aside class="notes">
                (~2 min) The money slide of the section: this is a real
                reconstruction from the checkpoint we trained for this talk, not
                an illustration. The model saw zeros in the shaded windows and
                drew the green curve. Look at what it got right, the line
                positions, and broadly the depths. It cannot have interpolated
                them; it had to infer "this is a star with these parameters, so
                these lines go here, this deep". That inference is the
                chemistry, and it is why the latent is a chemical-tagging
                feature. Then the crucial detail: the loss is computed ONLY on
                hidden pixels. Score the visible ones too and the model wins by
                copying its input. MSE 0.144 in standardised flux units, not
                perfect, and it should not be: a model that reconstructed the
                noise would be memorising, not generalising.
            </aside>
        </section>

        <!-- 47c-b5 · Step 4 — the latent -->
        <section>
            <div class="eyebrow">
                Self-supervised &middot; masked AE &middot; step 4 of 4
            </div>
            <h2>Throw the decoder away, keep z</h2>
            <div
                class="fig-split"
                style="
                    --cols: 1.15fr 1fr;
                    margin-top: 0.2em;
                    align-items: center;
                "
            >
                <div class="figure" style="margin: 0">
                    <img
                        :src="asset('mae_step4_latent.png')"
                        alt="PCA-2D of the 256-d latent for 378 real member stars: five clusters land in distinct groups"
                        style="width: 100%; height: auto; display: block"
                    />
                </div>
                <div>
                    <ul class="checklist small">
                        <li>
                            Training done,
                            <strong>discard the decoder</strong>
                        </li>
                        <li>
                            Every star &rarr; <code>embed(x)</code> &rarr;
                            <strong>256 numbers</strong>
                        </li>
                        <li>
                            No masking at inference: the full spectrum goes in
                        </li>
                        <li>Those 256 numbers replace the 16 abundances</li>
                        <li>
                            Cluster them with <em>any</em> tool from this
                            lecture
                        </li>
                    </ul>
                    <p class="small muted" style="margin-top: 0.4em">
                        378 real member stars, five clusters, plain PCA-2D of
                        the latent.
                        <strong
                            >The model was never told any of these
                            labels.</strong
                        >
                    </p>
                </div>
            </div>
            <aside class="notes">
                (~1.5 min) Close the loop back to the lecture's spine. Training
                is a means, not the end: we keep the encoder, drop the decoder,
                and stop masking, at inference the whole spectrum goes in and
                256 numbers come out. Those numbers are a drop-in replacement
                for the 16 abundances, which means every algorithm from the last
                ninety minutes still applies: K-means, HDBSCAN*, UMAP, EVoC,
                unchanged, just on a different feature vector. The figure is the
                reward: five real clusters, plain linear PCA of the latent, and
                they land in distinct regions. Nothing in training knew that M 3
                and M 67 exist. Now the fair question, is this actually better
                than the abundances, and better than PCA on the same pixels?
                That is the benchmark, next.
            </aside>
        </section>

        <!-- 47c-c · Why the objective matters -->
        <section>
            <div class="eyebrow">Self-supervision &middot; the objective</div>
            <h2>Reconstruction spends the budget on chemistry</h2>
            <div class="cols" style="--n: 2; margin-top: 0.4em">
                <div class="panel">
                    <h3>PCA, maximise variance</h3>
                    <ul class="dotlist small">
                        <li>
                            objective: keep the directions of
                            <strong>largest pixel variance</strong>
                        </li>
                        <li>budget spent on the continuum + temperature</li>
                        <li>
                            <strong>linear</strong>, one global linear map for
                            all stars
                        </li>
                        <li>
                            the element lines are low-variance &rarr; discarded
                        </li>
                    </ul>
                </div>
                <div class="panel flip">
                    <h3>Masked AE, maximise predictability</h3>
                    <ul class="dotlist small">
                        <li>
                            objective:
                            <strong>reconstruct hidden blocks</strong> from the
                            latent
                        </li>
                        <li>
                            budget spent on the local line physics (the only way
                            to fill a block)
                        </li>
                        <li>
                            <strong>nonlinear</strong>, a learned manifold per
                            stellar regime
                        </li>
                        <li>the lines are the only clue &rarr; preserved</li>
                    </ul>
                </div>
            </div>
            <p class="small muted center" style="margin-top: 0.4em">
                <strong>Variance &ne; information.</strong> The continuum has
                the variance; the chemistry has the physics. PCA keeps the
                first; the masked AE is forced to learn the second. (And on the
                re-run the AE still clears the abundances (UMAP 0.76 vs 0.58)
                even though PCA is a surprisingly strong linear baseline.)
            </p>
            <aside class="notes">
                (~2 min) The conceptual core. Two different objectives answer
                two different questions. PCA answers "which directions carry the
                most pixel variance?", the continuum, the temperature, the
                reddening. The masked AE answers "which latent lets me predict
                the pixels I hid?", and the only way to predict a hidden Fe I
                block is to know the iron, the temperature, and the blends.
                Variance is cheap and global; predictability is expensive and
                local. This is why the same 256-d latent is qualitatively
                different: one is a variance summary, the other is a predictive
                physical model. And the nonlinearity matters: the
                spectrum&rarr;abundance map is not a linear subspace, so a
                linear projector cannot align its axes with the chemistry. One
                honesty note for the Q&A: on the uniform re-run PCA is
                empirically competitive (0.73/0.75/0.74 vs our 0.69/0.76/0.74),
                so the claim is not "we crush PCA"; it is "we beat the
                abundances, and the objective is the right one".
            </aside>
        </section>

        <!-- 47c-c2 · The five clusters -->
        <section>
            <div class="eyebrow">Head-to-head &middot; the subjects</div>
            <h2>Five real clusters, one hard question</h2>
            <div class="cols" style="--n: 5; margin-top: 0.35em; gap: 0.5em">
                <div
                    class="panel"
                    v-for="c in [
                        {
                            img: 'sky_berkeley66.png',
                            name: 'Berkeley 66',
                            kind: 'open',
                            dist: '5.3 kpc',
                            n: '20',
                        },
                        {
                            img: 'sky_ic166.png',
                            name: 'IC 166',
                            kind: 'open',
                            dist: '4.9 kpc',
                            n: '13',
                        },
                        {
                            img: 'sky_m3.png',
                            name: 'M 3',
                            kind: 'globular',
                            dist: '10.2 kpc',
                            n: '142',
                        },
                        {
                            img: 'sky_m67.png',
                            name: 'M 67',
                            kind: 'open',
                            dist: '0.86 kpc',
                            n: '271',
                        },
                        {
                            img: 'sky_ngc188.png',
                            name: 'NGC 188',
                            kind: 'open',
                            dist: '1.9 kpc',
                            n: '27',
                        },
                    ]"
                    :key="c.name"
                >
                    <div
                        class="figure"
                        style="aspect-ratio: 1 / 1; margin: 0 0 0.25em 0"
                    >
                        <img
                            :src="asset(c.img)"
                            :alt="'DSS2 sky image of ' + c.name"
                            style="
                                width: 100%;
                                height: 100%;
                                object-fit: cover;
                                border-radius: 4px;
                            "
                        />
                    </div>
                    <h3 style="margin: 0">{{ c.name }}</h3>
                    <p
                        class="small muted"
                        style="margin: 0.05em 0 0 0; font-size: 0.42em"
                    >
                        {{ c.kind }} &middot; {{ c.dist }} &middot;
                        {{ c.n }} known members
                    </p>
                </div>
            </div>
            <p class="small muted center" style="margin-top: 0.35em">
                DSS2 sky images (SkyView, NASA GSFC). One globular (M 3) and
                four open clusters; the globular's presence is exactly what the
                control on the next slide is for.
            </p>
            <aside class="notes">
                (~45 s) Meet the subjects before the numbers. These are the five
                clusters the head-to-head table scores: four open clusters plus
                one globular, M 3. Real images, real distances, real member
                counts from the literature catalogues. Keep one fact in mind: M
                3 is a globular at 10 kpc, completely different metallicity
                regime, so "we separate the five clusters" could be one easy
                split. We will test exactly that in a moment.
            </aside>
        </section>

        <!-- 47c-c3 · The kinematic referee -->
        <section>
            <div class="eyebrow">Head-to-head &middot; the labels</div>
            <h2>The referee is kinematics, not chemistry</h2>
            <div
                class="fig-split"
                style="
                    --cols: 1.25fr 1fr;
                    margin-top: 0.35em;
                    align-items: center;
                "
            >
                <div class="figure" style="aspect-ratio: 1675 / 934">
                    <img
                        :src="asset('proper_motions.png')"
                        alt="Proper motions: each cluster's members form a tight clump in the Gaia proper-motion plane, the field stars are spread"
                        style="width: 100%; height: 100%; object-fit: contain"
                    />
                </div>
                <div>
                    <ul class="checklist small">
                        <li>
                            <strong
                                >Labels come from Gaia proper motions</strong
                            >
                            + parallax + radial velocity, who moves together
                            belongs together.
                        </li>
                        <li>
                            <strong>No chemistry in the ground truth.</strong>
                            The benchmark asks whether the latent can recover a
                            label it never saw, and the label was not made with
                            abundances.
                        </li>
                        <li>
                            Members are the magenta clumps; the field is
                            everything else in the cluster's patch of sky.
                        </li>
                    </ul>
                </div>
            </div>
            <p class="small muted center" style="margin-top: 0.35em">
                A clean experiment:
                <strong
                    >unsupervised features → clusters → kinematic
                    referee</strong
                >. Nothing in the loop touches element ratios.
            </p>
            <aside class="notes">
                (~1 min) The benchmark's labels are kinematic: Gaia proper
                motions (the panels above), parallax, and APOGEE radial
                velocities. Cluster members form a tight clump in the
                proper-motion plane; the magenta dots; the field is the grey
                haze. This matters because the whole claim is "spectra know
                chemistry, and chemistry knows clusters." If our labels had been
                built from abundances, the benchmark would be circular. They are
                not: the referee and the tested features are physically
                independent.
            </aside>
        </section>

        <!-- 47c-d · The head-to-head -->
        <section>
            <div class="eyebrow">Head-to-head &middot; DR19 uniform re-run</div>
            <h2>
                The masked latent beats the abundances, PCA gets a fair shot too
            </h2>
            <div class="panel">
                <h3>
                    Cluster-only homogeneity, same 982 stars, same 25 clusters,
                    mean &plusmn; sd over 7 seeds
                </h3>
                <table style="font-size: 0.5em; margin-top: 0.25em">
                    <thead>
                        <tr>
                            <th>features (all unsupervised)</th>
                            <th>t-SNE</th>
                            <th>UMAP</th>
                            <th>EVoC</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>abundances (16-d)</td>
                            <td>0.56 &plusmn; 0.00</td>
                            <td>0.58 &plusmn; 0.01</td>
                            <td>0.42 &plusmn; 0.04</td>
                        </tr>
                        <tr>
                            <td>PCA 64-d (linear)</td>
                            <td>0.74 &plusmn; 0.00</td>
                            <td>0.75 &plusmn; 0.01</td>
                            <td>0.73 &plusmn; 0.01</td>
                        </tr>
                        <tr>
                            <td>PCA 256-d (linear)</td>
                            <td>0.69 &plusmn; 0.00</td>
                            <td>0.69 &plusmn; 0.01</td>
                            <td>0.65 &plusmn; 0.02</td>
                        </tr>
                        <tr>
                            <td><strong>masked AE 256-d</strong></td>
                            <td><strong>0.74 &plusmn; 0.00</strong></td>
                            <td><strong>0.76 &plusmn; 0.01</strong></td>
                            <td><strong>0.69 &plusmn; 0.02</strong></td>
                        </tr>
                    </tbody>
                </table>
                <p
                    class="small muted"
                    style="margin-top: 0.3em; font-size: 0.42em"
                >
                    All four arms now embed the
                    <strong>same raw DR19 mwmStar spectra</strong> through the
                    re-run (see the confound slide). Every row is the
                    intersection of the arms' star lists.
                </p>
            </div>
            <p class="small muted center" style="margin-top: 0.4em">
                <strong>UMAP 0.76 vs 0.58</strong> for ASPCAP abundances on the
                identical sample, a real, seed-stable gap. But
                <strong>PCA 64-d is now competitive</strong> (0.75 UMAP): the
                old &ldquo;3&times; the linear baseline&rdquo; was a
                product-mismatch artefact, and we say so.
            </p>
            <aside class="notes">
                (~1.5 min) The honest head-to-head, re-run on uniform raw DR19
                mwmStar after we caught the product mismatch. Every row is the
                same 982 stars across 25 clusters, same clustering, error bars
                over 7 seeds. Two things survive, and one does not. Survives:
                the masked AE beats the ASPCAP abundances clearly, UMAP 0.76 vs
                0.58, EVoC 0.69 vs 0.42, and that gap is far outside the seed
                spread, so "a model that never saw an element ratio separates
                clusters better than the element ratios do" still holds. Does
                not survive: the claim that we are 3x better than PCA. On the
                full sample PCA-64d is 0.73/0.75/0.74 against our 0.69/0.76/0.74
                essentially tied. The earlier "3x PCA" came from the product
                mismatch inflating our arm and the degenerate PCA arm
                collapsing. Be upfront: the self-supervised win over a linear
                baseline is modest; the real win is over abundances, and over
                PCA only on EVoC.
            </aside>
        </section>

        <!-- 47c-d1 · See it with your eyes -->
        <section class="densest">
            <div class="eyebrow">Head-to-head &middot; look at it</div>
            <h2>Same stars, two spaces, which one knows the clusters?</h2>
            <div
                class="figure"
                style="
                    aspect-ratio: 1456 / 765;
                    width: 82%;
                    max-height: 56vh;
                    margin: 0.25em auto 0;
                "
            >
                <img
                    :src="asset('headtohead_pca.png')"
                    alt="Plain linear 2-D view of the same stars: in the masked AE latent five of the clusters are compact isolated islands; in the abundances they are loose and overlap"
                    style="width: 100%; height: 100%; object-fit: contain"
                />
            </div>
            <p class="small muted center" style="margin-top: 0.3em">
                No t-SNE trickery, a
                <strong>plain linear 2-D projection</strong> of each space.
                <strong>Left:</strong> the masked AE latent;
                <strong>right:</strong> the 16 ASPCAP abundances on the same 378
                stars (five of the 25 clusters). The silhouette scores in the
                panels (0.62 vs &minus;0.07) count what your eyes see.
            </p>
            <aside class="notes">
                (~1 min) The table, rendered as a picture, and deliberately
                without any nonlinear projection, because t-SNE would force
                separation in both panels and lie to you. This is PCA-2D of each
                space on 378 stars across five of the 25 clusters: on the left
                the clusters are compact islands; on the right they are loose,
                overlapping clouds, the abundance silhouette goes slightly
                negative, meaning members are closer to other clusters than to
                their own. The silhouette score (0.62 vs -0.07) quantifies
                exactly the visual difference. Same stars, same clusters, no
                distortion; this is the 0.76 vs 0.58 in the table, made visible.
            </aside>
        </section>

        <!-- 47c-d2 · Is it just the globular? -->
        <section>
            <div class="eyebrow">Head-to-head &middot; the control</div>
            <h2>Not just &ldquo;globular vs open&rdquo;</h2>
            <p class="small">
                M 3 is the largest single cluster in the sample. A globular sits
                at a completely different metallicity, so &ldquo;we separate
                clusters&rdquo; could just mean &ldquo;we spotted the
                globular.&rdquo; So drop it and re-score the
                <strong>24 clusters that remain</strong> (884 stars), the
                genuinely hard case.
            </p>
            <div class="panel" style="margin-top: 0.4em">
                <h3>
                    24 clusters, M 3 removed (the same 884 stars for both rows,
                    7 seeds)
                </h3>
                <table style="font-size: 0.5em; margin-top: 0.25em">
                    <thead>
                        <tr>
                            <th>features</th>
                            <th>t-SNE</th>
                            <th>UMAP</th>
                            <th>EVoC</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>abundances (16-d)</td>
                            <td>0.20 &plusmn; 0.00*</td>
                            <td>0.58 &plusmn; 0.01</td>
                            <td>0.50 &plusmn; 0.03</td>
                        </tr>
                        <tr>
                            <td><strong>masked AE 256-d</strong></td>
                            <td><strong>0.71 &plusmn; 0.01</strong></td>
                            <td><strong>0.74 &plusmn; 0.01</strong></td>
                            <td><strong>0.66 &plusmn; 0.02</strong></td>
                        </tr>
                    </tbody>
                </table>
                <p
                    class="small muted"
                    style="font-size: 0.42em; margin: 0.25em 0 0"
                >
                    * collapsed: on the abundances, t-SNE → HDBSCAN* finds only
                    2 clusters, with 75% of the stars in one, a floor, not a
                    comparison. Compare on UMAP and EVoC.
                </p>
            </div>
            <p class="small muted center" style="margin-top: 0.4em">
                The masked AE holds at 0.71/0.74/0.66, barely below the full
                sample, and still clears the abundances (0.74 vs 0.58 UMAP). The
                globular was never doing the work.
            </p>
            <aside class="notes">
                (~1 min) The control that matters. Remove the globular, keep 24
                clusters across the full metallicity range, and the masked
                latent barely moves (0.71/0.74/0.66 vs 0.74/0.76/0.69 on the
                full sample) while still clearing the abundances (0.74 vs 0.58
                UMAP). Both rows are scored on exactly the same 884 stars. Do
                not quote the abundance t-SNE number as a win for the AE: that
                partition collapsed to two clusters, so 0.20 says the clusterer
                failed, not that the features are three times worse. The honest
                version: the earlier "3.7x on 31 open-cluster stars" was a small
                sample inflated by the product mismatch; the re-run says the
                AE's edge over abundances is ~1.3x but it survives dropping the
                globular and is seed-stable. If someone asks "isn't this just
                finding M 3?", the answer is no, the number barely changes
                without it.
            </aside>
        </section>

        <!-- 47c-d3 · The full experiment -->
        <section class="densest">
            <div class="eyebrow">Head-to-head &middot; the full experiment</div>
            <h2>Not five clusters, twenty-five</h2>
            <div
                class="fig-split"
                style="
                    --cols: 1.35fr 1fr;
                    margin-top: 0.35em;
                    align-items: center;
                "
            >
                <div
                    class="figure"
                    style="aspect-ratio: 1530 / 1530; max-height: 52vh"
                >
                    <img
                        :src="asset('sky_montage_25.png')"
                        alt="Montage of all 25 clusters: 18 open clusters and 7 globulars, DSS2 images"
                        style="width: 100%; height: 100%; object-fit: contain"
                    />
                </div>
                <div>
                    <p class="small">
                        The whole catalogue:
                        <strong>18 open clusters + 7 globulars</strong>, every
                        one with DSS2 imaging and Gaia kinematics. Cluster-only
                        homogeneity on all members (uniform re-run):
                    </p>
                    <div class="panel" style="margin-top: 0.35em">
                        <table style="font-size: 0.5em; margin-top: 0.2em">
                            <thead>
                                <tr>
                                    <th>features (25 clusters, 982 stars)</th>
                                    <th>t-SNE</th>
                                    <th>UMAP</th>
                                    <th>EVoC</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td>abundances (16-d)</td>
                                    <td>0.56</td>
                                    <td>0.58</td>
                                    <td>0.42</td>
                                </tr>
                                <tr>
                                    <td><strong>masked AE 256-d</strong></td>
                                    <td><strong>0.74</strong></td>
                                    <td><strong>0.76</strong></td>
                                    <td><strong>0.69</strong></td>
                                </tr>
                                <tr>
                                    <td>kinematics only (4-d)</td>
                                    <td>0.94</td>
                                    <td>0.94</td>
                                    <td>0.88</td>
                                </tr>
                            </tbody>
                        </table>
                        <p
                            class="small muted"
                            style="font-size: 0.38em; margin-top: 0.3em"
                        >
                            Field retrieval, one shared star list (30 107 stars,
                            978 kinematic members, chance precision 3%):
                            abundances t-SNE recall 0.20 / precision 0.24, UMAP
                            1.00 / 0.001 (a blob), EVoC 0.56 / 0.003; the
                            spectral arms (PCA 64-d, masked AE) reach recall
                            0.42 / precision 0.62. Spectra double the recall and
                            lift precision 2.6×, real, but still not a fishing
                            net.
                        </p>
                    </div>
                </div>
            </div>
            <aside class="notes">
                (~1 min) The scope claim: this is not a five-cluster anecdote,
                it is the full 25-cluster benchmark; every cluster with real sky
                imaging and kinematic membership. Cluster-only homogeneity holds
                at scale: 0.76 vs 0.58 (UMAP) on 982 stars. And now that the
                re-run is uniform, field retrieval is measurable: t-SNE recall
                0.21 / precision 0.22 against ~3% chance, honest, and modest.
                Two things to be transparent about: we are still 0.76 vs 0.94
                behind kinematics (physics says we should be), and the abundance
                arm here is the same pipeline on the same stars, so this row IS
                comparable.
            </aside>
        </section>

        <!-- 47c-e · Takeaway -->
        <section>
            <div class="eyebrow">Self-supervision &middot; the takeaway</div>
            <h2>Dimensionality reduction, done by physics</h2>
            <p class="small">
                Both PCA and the masked AE compress the same 8575 pixels into a
                compact feature. The difference is the
                <strong>objective</strong>:
            </p>
            <div class="cols" style="--n: 2; margin-top: 0.4em">
                <div class="panel">
                    <h3>PCA</h3>
                    <p class="small">
                        &ldquo;keep the directions of largest
                        variance.&rdquo;<br />A linear summary of the pixels.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Masked AE</h3>
                    <p class="small">
                        &ldquo;predict the parts I hid from you.&rdquo;<br />A
                        latent that must model the line physics.
                    </p>
                </div>
            </div>
            <p class="small muted center" style="margin-top: 0.4em">
                <strong>The spectrum is its own label.</strong> No catalogue, no
                element ratios, no circularity, and the resulting latent is the
                best chemical-tagging feature we measured.
            </p>
            <aside class="notes">
                (~1 min) Close the arc. "Dimensionality reduction" is
                underspecified; it is the objective that decides what the
                reduced space keeps. PCA keeps variance; the masked AE keeps
                predictability. Because the element ratios are a low-variance,
                high-information signal, the variance objective throws them away
                while the reconstruction objective cannot avoid them. The
                spectrum is its own label: no ASPCAP catalogue, no circularity,
                and the latent beats everything we measured. This reframes the
                workshop's own method; the RNN's real value was never the
                regression head, it was the learned compression.
            </aside>
        </section>

        <!-- 47d · Two surveys -->
        <section>
            <div class="eyebrow">
                Two surveys &middot; giants + main sequence
            </div>
            <h2>Distance from the clump, age from the turnoff</h2>
            <p class="small">
                GALAH (DEC &#8818; +25&deg;) sees the
                <strong>main sequence + turnoff</strong>; APOGEE sees the
                <strong>giants</strong>. Combine them, 6 clusters, 14 common
                abundances, Gaia parallax/PM/RV, and the cluster parameters
                finally separate: the red clump pins the distance, the main
                sequence pins the age.
            </p>
            <div class="cols" style="--n: 2; margin-top: 0.4em">
                <div class="panel">
                    <h3>Red-clump distance (J&minus;K)</h3>
                    <table style="font-size: 0.5em; margin-top: 0.25em">
                        <thead>
                            <tr>
                                <th>cluster</th>
                                <th>dm clump</th>
                                <th>dm lit</th>
                                <th>resid</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>NGC 6819</td>
                                <td>11.92</td>
                                <td>11.90</td>
                                <td><strong>+0.01</strong></td>
                            </tr>
                            <tr>
                                <td>NGC 2243</td>
                                <td>13.08</td>
                                <td>13.25</td>
                                <td><strong>&minus;0.16</strong></td>
                            </tr>
                            <tr>
                                <td>NGC 7789</td>
                                <td>11.68</td>
                                <td>11.27</td>
                                <td>+0.41</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
                <div class="panel flip">
                    <h3>NGC 2243, age over 6 seeds</h3>
                    <table style="font-size: 0.5em; margin-top: 0.25em">
                        <thead>
                            <tr>
                                <th></th>
                                <th>stars</th>
                                <th>age, Gyr (lit 1.08)</th>
                                <th>closer in</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>APOGEE giants</td>
                                <td>38</td>
                                <td>1.02 – 1.59</td>
                                <td>3 of 6</td>
                            </tr>
                            <tr>
                                <td>MS + giants</td>
                                <td>256</td>
                                <td>0.64 – 0.98</td>
                                <td>3 of 6</td>
                            </tr>
                        </tbody>
                    </table>
                    <p
                        class="small muted"
                        style="margin-top: 0.25em; font-size: 0.42em"
                    >
                        Clump distance used as a &plusmn;0.2 mag prior. Adding
                        the main sequence moves the age down, but not reliably
                        closer to the literature.
                    </p>
                </div>
            </div>
            <p class="small muted center" style="margin-top: 0.4em">
                Chemical tagging is a <strong>two-survey problem</strong>:
                Kos&rsquo;s GALAH saw the turnoff (age), Garcia-Dias&rsquo;s
                APOGEE saw the giants (distance). The clump distance already
                works; the joint age fit does not yet.
            </p>
            <aside class="notes">
                (~2 min) The synthesis that closes the loop on the two opening
                papers. GALAH reaches only DEC &lt; +25 deg, so we added two
                southern intermediate-age clusters (NGC 2243, Collinder 261)
                where both surveys overlap. The red-clump distance is the median
                K of the J-K clump slice (member giants, metallicity-corrected
                M_K): NGC 6819 lands within 0.01 mag, NGC 2243 within 0.16 mag.
                Then the sweet spot: the APOGEE member giants locate the
                cluster, the GALAH main sequence is selected around that
                kinematic centroid with a relaxed parallax (at 4.7 kpc the
                parallax error rivals the parallax). The clump numbers are
                deterministic and reproduce exactly on DR19. The joint isochrone
                fit does not hold up. An earlier single run gave 0.98 Gyr
                against 1.59 for the giants alone, but across six seeds the
                giants-only fit ranges 1.02–1.59 Gyr and the combined fit
                0.64–0.98, and each is closer to the literature 1.08 in three of
                the six. The age posterior (sd 0.65–1.04 dex) is close to the
                width of its flat prior (1.0) either way. So say it plainly: the
                idea is right (giants carry distance, the main-sequence turnoff
                carries age), the distance half works, and the age half needs a
                better likelihood before it shows anything.
            </aside>
        </section>

        <!-- 47e · Literature — the best-case sample, published -->
        <section class="dense">
            <div class="eyebrow">
                Positioning &middot; the best-case sample, published
            </div>
            <h2>The same experiment, published</h2>
            <div class="cols" style="--n: 2; margin-top: 0.4em">
                <div class="panel">
                    <h3>Casamiquela et al. 2021; what they had</h3>
                    <ul
                        class="small"
                        style="margin: 0.25em 0 0 1em; line-height: 1.4"
                    >
                        <li>
                            175 red-clump members of
                            <strong>31 thin-disc open clusters</strong>, and no
                            field stars to confuse the question.
                        </li>
                        <li>
                            16 elements measured
                            <strong>differentially</strong>: intra-cluster
                            spread typically <strong>0.03 dex</strong>,
                            individual errors usually &le; 0.05 dex.
                        </li>
                        <li>
                            HDBSCAN straight on the abundances, free parameters
                            fine-tuned to
                            <strong
                                >maximise the number of clusters
                                recovered</strong
                            >
                            a sample they call their own &ldquo;best-case
                            scenario&rdquo;.
                        </li>
                    </ul>
                </div>
                <div class="panel flip">
                    <h3>What they got</h3>
                    <ul
                        class="small"
                        style="margin: 0.25em 0 0 1em; line-height: 1.4"
                    >
                        <li>
                            <strong
                                >h 0.49 &middot; c 0.63 &middot; V 0.55</strong
                            >, and <strong>9 of 31</strong> clusters recovered
                            at 40 %; at 70 %, one: <strong>NGC 2682</strong>.
                        </li>
                        <li>
                            More than
                            <strong
                                >70 % of the groups found are
                                statistical</strong
                            >
                            they mix stars from different clusters. Four of the
                            nine recovered groups are 100 % pure, and most of
                            those hold only two of the cluster's stars.
                        </li>
                        <li>
                            Put the field back in (APOGEE DR16 red clump, 16 193
                            stars, 18 astroNN species &rarr; 272 groups): only
                            <strong>NGC 188</strong> clears 40 %, nothing clears
                            70 %, and cutting the element list does not help.
                        </li>
                    </ul>
                </div>
            </div>
            <p class="small muted" style="margin: 0.4em 0 0 170px">
                Why it cannot go better, in their words: the clusters' chemical
                signatures
                <em>overlap widely</em>, across only &minus;0.2 to 0.1 dex in
                [X/H]; the disc's birth gas was <strong>well mixed</strong>,
                with a scatter of 0.02–0.03 dex (Kreckel et al. 2020) that is
                the size of the measurement precision itself.
            </p>
            <aside class="notes">
                (~1.5 min) The paper to beat, and the closest published
                measurement of the question this talk asks. Their sample is the
                best case spectroscopy can buy: differential abundances, 16
                elements, per-element precision three-fifths of ours (0.03 dex
                against our ~0.05), no field stars, and clustering parameters
                tuned for recovery rather than for honesty. It still recovers 9
                of 31 clusters, and their own field arm echoes our field slide:
                272 groups, one cluster above 40 %. State the mechanism plainly,
                because it answers "why is this hard": the clusters formed from
                a well-mixed disc, so the chemical difference between two
                clusters is the size of the error bar. Then flag the two
                differences we cannot control and do not hide them; they have ~6
                stars per cluster from high-resolution spectroscopy; we have ~39
                per cluster from survey ASPCAP abundances.
            </aside>
        </section>

        <!-- 47e2 · Literature — their metric, our stars -->
        <section class="densest">
            <div class="eyebrow">
                Positioning &middot; their metric, our stars
            </div>
            <h2>Where we land against it</h2>
            <table style="font-size: 0.5em; margin-top: 0.3em">
                <thead>
                    <tr>
                        <th>run</th>
                        <th>stars / clusters</th>
                        <th>V</th>
                        <th>chance V</th>
                        <th>RF<sub>40</sub></th>
                        <th>RF<sub>70</sub></th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>Casamiquela+21, published</td>
                        <td>175 / 31</td>
                        <td>0.55</td>
                        <td>0.44–0.51</td>
                        <td><strong>0.29</strong></td>
                        <td>0.03</td>
                    </tr>
                    <tr>
                        <td>their clustering step, our stars (raw dex)</td>
                        <td>982 / 25</td>
                        <td>0.45</td>
                        <td>0.29</td>
                        <td>0.08</td>
                        <td>0.04</td>
                    </tr>
                    <tr>
                        <td>ours: t-SNE on 16 abundances</td>
                        <td>982 / 25</td>
                        <td>0.58</td>
                        <td>0.09</td>
                        <td>0.32</td>
                        <td>0.08</td>
                    </tr>
                    <tr>
                        <td>ours: UMAP on 16 abundances</td>
                        <td>982 / 25</td>
                        <td><strong>0.58</strong></td>
                        <td>0.11</td>
                        <td>0.33</td>
                        <td>0.01</td>
                    </tr>
                    <tr>
                        <td>ours: UMAP, open clusters only</td>
                        <td>665 / 18</td>
                        <td>0.49</td>
                        <td>0.13</td>
                        <td><strong>0.40</strong></td>
                        <td>0.04</td>
                    </tr>
                </tbody>
            </table>
            <ul class="small" style="margin: 0.45em 0 0 1em; line-height: 1.45">
                <li>
                    <strong
                        >V-measure does not transfer between samples.</strong
                    >
                    Shuffle their labels into groups of the observed sizes and a
                    random partition already scores V = 0.44–0.51, so their
                    published 0.55 sits 0.04–0.11 above their own chance level;
                    our 0.58 sits 0.45–0.49 above ours. Recovery fraction is the
                    comparison that does carry; its chance level is &le; 0.04 in
                    both samples.
                </li>
                <li>
                    <strong
                        >Their pipeline recovers less on our stars than on
                        theirs</strong
                    >
                    (0.08 against 0.29); ours lands in their range (0.32–0.40)
                    on survey-quality abundances.
                </li>
                <li>
                    <strong
                        >The survivors are the edge cases in both
                        papers.</strong
                    >
                    Theirs: NGC 2420, the most metal-poor cluster in their
                    sample, plus NGC 6705 and NGC 2682. Ours: NGC 2420 again,
                    recovered in 7 of 7 seeds by UMAP on the open clusters.
                </li>
            </ul>
            <p class="small muted" style="margin: 0.4em 0 0 170px">
                Not a head-to-head: different samples and different abundance
                quality. One published pipeline and its metric triple, not a
                survey of the field, and we do not re-run the Spina et al.
                graph-attention autoencoder (that comparison is on the
                &ldquo;two papers&rdquo; slide).
            </p>
            <aside class="notes">
                (~2 min) The point of this slide is that two V-measures from two
                samples are not the same number, and we say so before anyone
                else does. V depends on how many groups you cut and how big they
                are, so a random partition of their 175 stars already scores
                0.44–0.51: their published 0.55 is 0.04–0.11 above their floor,
                while ours is 0.45–0.49 above ours. Different quantity,
                different partition, different sample. What does carry across is
                the recovery fraction: 9 of 31 for them; 2 of 25 when we run
                their pipeline on our stars; 8 of 25 for our t-SNE and UMAP
                arms; 7 of 18 open clusters for the open-only run, which is the
                subset closest to their sample. The honest summary is "we match
                their best case on coarser abundances", not "we beat them". Then
                the caveats, out loud: ~6 stars per cluster against our ~39,
                differential abundances against survey ASPCAP, one published
                pipeline rather than a survey of the field. If someone asks why
                their pipeline drops to 0.08 on our stars, the two candidate
                reasons are the abundance precision and the mixed evolutionary
                stages in our sample, and we cannot separate them yet.
            </aside>
        </section>

        <!-- 48 · Lessons -->
        <section class="dense">
            <div class="eyebrow">What the numbers teach</div>
            <h2>Eight lessons from the benchmark</h2>
            <div class="slide-body">
                <ol class="contribs tight small">
                    <li>
                        <strong>All-sky collapses.</strong> One embedding of
                        every clean star buries each cluster in the field, cut a
                        sky region around it, as Kos et al. did (40° around the
                        Pleiades).
                    </li>
                    <li>
                        <strong>Precision is the hard part.</strong> Even region
                        by region, macro precision is ≈ 0.15 (t-SNE) and
                        0.02–0.06 (UMAP, EVoC): the field is chemically similar.
                        Kinematics confirm what chemistry only suggests.
                    </li>
                    <li>
                        <strong>Recall ≈ 1.0 is usually a blob.</strong>
                        HDBSCAN* swallows the dense field into one cluster, so
                        UMAP "recovers" everything at precision ≈ 0. kNN purity
                        is the honest, parameter-free score.
                    </li>
                    <li>
                        <strong
                            >Globulars tag cleanly, open clusters don't.</strong
                        >
                        M 3 / M 5 / M 15 reach kNN purity ≈ 0.4–0.5;
                        solar-metallicity open clusters stay ≲ 0.2, the paper's
                        own caveat.
                    </li>
                    <li>
                        <strong
                            >t-SNE isolates, UMAP over-merges, EVoC is fast but
                            coarse.</strong
                        >
                        Three sets of assumptions, three answers.
                    </li>
                    <li>
                        <strong>Levers are fragile.</strong> On DR17, weighting
                        each element by 1/σ <em>after</em> standardising lifted
                        M 67's t-SNE recall 0.09 → 0.44. On DR19 the same switch
                        moves recall 0.10 → 0.96 but precision 0.14 → 0.05:
                        t-SNE's cluster grows to swallow the field. Always read
                        recall and precision together.
                    </li>
                    <li>
                        <strong>Spectra beat abundances, but not PCA.</strong>
                        A self-supervised masked autoencoder on the raw spectrum
                        separates the 25 clusters better than the 16 ASPCAP
                        abundances (0.76 vs 0.58, UMAP, the same 982 stars) and
                        a plain PCA of the same spectra does about as well
                        (0.75).
                    </li>
                    <li>
                        <strong>Two surveys, two parameters.</strong> APOGEE
                        sees giants (red clump → distance, NGC 6819 within 0.01
                        mag); GALAH sees the main sequence (turnoff → age). The
                        joint age fit is not stable yet: it moves more between
                        seeds than between samples.
                    </li>
                </ol>
                <p class="small muted center" style="margin-top: 0.5em">
                    And one quieter lesson: row-normalisation is what breaks the
                    blob, precision rises three- to fivefold and recall falls.
                    No setting hands you both.
                </p>
            </div>
            <aside class="notes">
                (~3 min) The teaching payload; each lesson maps back to
                something in the tour. (1) is the region detail from the t-SNE
                slides. (2) is the paper's own "untaggable" caveat, 47 Tuc
                included. (3) is the validation slide's warning about scores
                that lie, and the reason we carry kNN purity at all. (4) is no
                free lunch in astronomical clothing: the right tool depends on
                the shape of the data, and metal-poor globulars simply are a
                different shape in C-space. (5) is the EVoC caveat about the
                metric. (6) is the element-weight lever, and a warning with it:
                the same switch that helped on DR17 inflates t-SNE's cluster on
                DR19 (one run each, M 67, 30°, 5000 field stars, seed 42; the
                full region sweep agrees that precision drops, 0.83 → 0.24), so
                a lever has to be checked on both recall and precision every
                time the data changes. (8) is the two-survey idea, stated as far
                as the numbers go. Close on the footnote: the precision/recall
                trade-off is a knob, not a bug; you choose which error you can
                live with.
            </aside>
        </section>

        <!-- 49 · Take-away -->
        <section>
            <div class="eyebrow">Take-away</div>
            <h2>No bad models, only mismatched ones</h2>
            <ul class="checklist medium">
                <li>
                    There are no bad models, only models applied outside their
                    assumptions
                </li>
                <li>Know your data first: shape, scale, density, metric</li>
                <li>
                    Statistics and visualisation matter as much as the model
                </li>
                <li>
                    Validate against independent evidence, and never quit
                    thinking
                </li>
            </ul>
            <p class="stat center" style="margin-top: 0.5em">
                Know the assumptions,<br />know your data
            </p>
            <aside class="notes">
                (~2 min) Land the through-line one more time, now earned by a
                concrete astronomy example. Walk back up the list on the slide:
                every clusterer you met today encodes assumptions about shape,
                scale, density and metric, and the failures you saw were never
                the algorithm being "bad"; they were the assumption being wrong
                for these data. The only way to know you chose right is to
                validate against something the model never saw. Here that was
                kinematics; in their own work it will be something else, but it
                has to exist.
            </aside>
        </section>

        <!-- 49a · Transition — part 2 → part 3. The 120-minute block starts
         here, before the assignment slide names the task and the hands-on slide
         gives the mechanics, so the ramp reads: statement → task → how. Kept
         deliberately short: 50 already carries the commands, the QR code and the
         doc pointers, so this slide only switches the room's mode and states the
         block's three-step shape. -->
        <section class="dense title-slide center">
            <div class="eyebrow">End of part 2 · 90 minutes done</div>
            <h1>Your turn.</h1>
            <p class="subtitle">
                Same data, same code, your own laptop, the last
                <strong>120 minutes</strong> are yours.
            </p>
            <div class="slide-body">
                <p class="small">
                    <strong>1 · Check it runs</strong>,
                    <code>./run.sh run --fast</code> (~2 min, no GPU)<br />
                    <strong>2 · The exercise</strong>, your cluster, from
                    <code>docs/cluster_assignment.md</code><br />
                    <strong>3 · Send it back</strong>, branch, push, open the
                    pull request
                </p>
                <p class="small muted" style="margin-top: 0.4em">
                    If your download is still finishing, tell me; we will sort
                    it out while the room gets started.
                </p>
            </div>
            <aside class="notes">
                (~1 min) Switch modes out loud here: from "here is what we
                found" to "here is what you do". The three steps are the whole
                block's shape (check, work, send back) and the next two slides
                give the task and the mechanics, so keep this one short and
                practical. Worth saying: how they get help (come to the front),
                and that the deliverable is a pull request, not a local
                notebook. Anyone whose download is still running should start at
                step 2 and check later.
            </aside>
        </section>

        <!-- 49a2 · Three ways to contribute -->
        <section class="dense">
            <div class="eyebrow">Your turn · pick a track</div>
            <h2>Three ways in, all of them real contributions</h2>
            <p class="small">
                The deliverable is a <strong>pull request</strong> either way,
                against
                <code>iaa-so-training/iaa-advanced-neural-networks-2026</code>,
                folder <code>day_4_clustering/</code>. Pick the one that suits
                how you like to work; they are worth the same.
            </p>
            <div class="cols" style="--n: 3; margin-top: 0.5em">
                <div class="panel">
                    <h3>1 · Beat the baseline</h3>
                    <p class="small">
                        Take your cluster, extend the pipeline, move recall,
                        precision or purity. The research track: open-ended,
                        and a null result honestly explained counts. Details on
                        the next slide.
                    </p>
                </div>
                <div class="panel">
                    <h3>2 · Improve the workbook</h3>
                    <p class="small">
                        The text you are reading is LaTeX in
                        <code>article/chapters/*.tex</code>, 17 chapters. Fix an
                        explanation that did not land, an error, a missing
                        citation into <code>references.bib</code>. If it
                        confused you, it will confuse the next reader.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>3 · Work the exercises</h3>
                    <p class="small">
                        <strong>56 exercises</strong> across 16 chapters, each
                        with a Jupyter deck and a worked solution module. Do
                        them, disagree with an answer, send a better one. The
                        guided track, and the slide after next.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.5em">
                Fork first, branch, then <em>Contribute → Open pull request</em>.
                Everything you need is in
                <code>day_4_clustering/CONTRIBUTING.md</code>.
            </p>
            <aside class="notes">
                (~1 min) Say plainly that the room is not one kind of person.
                Some want the open research problem, some would rather improve
                the writing, some learn by working problems with a solution to
                check against. All three land in the same repository through the
                same mechanism, a pull request, so all three teach the part
                this school actually cares about: contributing to someone
                else's codebase in public. Mention the fork requirement once
                here and then move on; they already forked it this morning for
                the download, and if they cloned the original instead,
                CONTRIBUTING has the section on repointing the remote rather
                than re-cloning five gigabytes. Then walk the three slides: the
                assignment next, then the workbook and exercises detail.
            </aside>
        </section>

        <!-- 49b · The assignment -->
        <section class="densest">
            <div class="eyebrow">The assignment · 25 clusters</div>
            <h2>Take a cluster, beat our baseline</h2>
            <p class="small">
                The results you just saw are a
                <strong>baseline, not a ceiling</strong>. Each of you gets one
                cluster from the region sweep,
                <code>docs/region_sweep_results.md</code> in the repo holds our
                current recall / precision / kNN-purity for every one. Your job:
                <strong
                    >extend the pipeline and get better numbers for your
                    cluster</strong
                >.
            </p>
            <div class="cols" style="--n: 3; margin-top: 0.45em">
                <div class="panel">
                    <h3>Your cluster</h3>
                    <p class="small">
                        18 open clusters (Pleiades, M 67, NGC 188, NGC 6819, NGC
                        2243, Collinder 261, …) and 7 globulars (M 3, M 5, M 13,
                        M 15, M 71, M 107, M 92). Pair up on the crowded fields.
                        The globulars already tag cleanly; your real job is the
                        open clusters that don't.
                    </p>
                </div>
                <div class="panel">
                    <h3>Levers we've already found</h3>
                    <p class="small">
                        Element-wise 1/σ weights (after standardising) the
                        single biggest gain. HDBSCAN
                        <code>min_cluster_size</code>, t-SNE
                        <code>perplexity</code>, SNR cut, dwarf-only, element
                        subset. Newer levers: the
                        <strong>RNN spectral latent</strong> (256-d), the
                        <strong>red-clump / isochrone</strong> fit, and the
                        <strong>GALAH main-sequence</strong> cross-match. All
                        logged in <code>docs/experiment_results.md</code>.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Deliverable</h3>
                    <p class="small">
                        One slide: your cluster, your best recall / precision /
                        purity, and the <em>one</em> change that moved them.
                        Beat the baseline, or explain honestly why your cluster
                        resists. Either answer is a result.
                    </p>
                </div>
            </div>
            <p class="stat center" style="margin-top: 0.5em">
                Not a bigger number,<br />knowing which knob moved why
            </p>
            <aside class="notes">
                (~2 min) This is the hinge between "here is what I did" and "now
                it is your turn". Hand out the cluster list, then make the
                framing explicit: the table in docs/region_sweep_results.md is
                the scoreboard, and the levers in docs/experiment_results.md are
                the map; on DR17 the biggest was element weights (0.09 → 0.44
                t-SNE recall on M 67), but on DR19 the same switch lowers
                precision in both of our re-runs, so even the known levers need
                re-checking on the current data. Their real job is the next
                lever, specific to their cluster. Stress that a null result is
                still a result: if a cluster refuses to tag cleanly, say why;
                that is the paper's own finding for 47 Tuc. Close by reminding
                them this is track 1 of three, and the other two are next.
            </aside>
        </section>

        <!-- 49c · Track 2, the workbook -->
        <section class="dense">
            <div class="eyebrow">Track 2 · the workbook</div>
            <h2>Send a PR to the text itself</h2>
            <p class="small">
                The companion workbook is LaTeX in the same repository, one file
                per chapter under <code>article/chapters/</code>, compiled into
                <code>workbook.pdf</code>. Prose is as reviewable as code, and a
                confusing paragraph is a defect.
            </p>
            <div class="cols compact" style="--n: 3; margin-top: 0.45em">
                <div class="panel">
                    <h3>What is worth a PR</h3>
                    <ul class="dotlist small" style="margin: 0.1em 0 0">
                        <li>A passage you had to read three times.</li>
                        <li>A wrong or missing number.</li>
                        <li>A claim cited to the wrong paper, or to none.</li>
                        <li>A figure whose caption does not say the point.</li>
                    </ul>
                </div>
                <div class="panel">
                    <h3>How to build it</h3>
                    <pre
                        class="small"
                        style="text-align: left; margin: 0.15em 0 0"
                    ><code>cd day_4_clustering/article
docker run --rm -v "$PWD:/w" -w /w \
  texlive/texlive latexmk -pdf workbook.tex</code></pre>
                    <p class="small" style="margin-top: 0.3em">
                        No LaTeX on your laptop. Not building at all is fine
                        too; the diff on the <code>.tex</code> is what gets
                        reviewed.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>The one hard rule</h3>
                    <p class="small">
                        Citations are <strong>keys</strong> into
                        <code>article/references.bib</code>, never a typed
                        author-year string. An unknown key fails the tests
                        instead of printing a reference nobody can follow. Add
                        the entry in the same PR.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.45em">
                17 chapters, 64 figures, one bibliography. Exercise statements
                live here too, so changing one is a change to the exercise.
            </p>
            <aside class="notes">
                (~1 min 30 s) This is the track for the people who would rather
                write than tune hyperparameters, and it is not the consolation
                prize: the workbook is what the next cohort reads, and a
                paragraph that lost you is a real defect in it. Point out that
                they are the ideal reviewers right now, today, because they have
                just met this material cold and will never be this unfamiliar
                with it again; by next week they will read past the confusing
                paragraph without noticing. Name the citation rule explicitly
                because it is the one thing that will bounce their PR: the bib
                is the single source of truth, cite() raises on a key that is
                not in it, and the test suite checks. The build command is a
                one-off container, the same trick as the workshop image: no
                TeX distribution on their laptop, and the figures are already
                committed so it compiles straight from a clean checkout
                (measured: 67 pages, no errors). Finally, warn them the
                exercise statements are in these same chapter files, so editing
                one means the exercise deck must be rebuilt, which is the next
                slide.
            </aside>
        </section>

        <!-- 49d · Track 3, the exercises -->
        <section class="denser">
            <div class="eyebrow">Track 3 · the exercises</div>
            <h2>56 exercises, 16 notebooks, one module each</h2>
            <p class="small">
                Every chapter has a Jupyter deck. Open it, work the problem in
                the empty cell, then reveal the worked answer. The notebooks
                <strong>compute nothing</strong>; each answer lives in a Python
                module you can open, read and re-run.
            </p>
            <div class="cols compact" style="--n: 2; margin-top: 0.4em">
                <div class="panel">
                    <h3>Open one</h3>
                    <pre
                        class="small"
                        style="text-align: left; margin: 0.15em 0 0"
                    ><code>cd day_4_clustering
mkdir -p data results notebooks
docker compose up      # localhost:9999</code></pre>
                    <p class="small" style="margin-top: 0.35em">
                        Then open
                        <code>notebooks/exercises/chapter_09_validation.ipynb</code>,
                        or <code>workbook_exercises.ipynb</code> for all 16.
                        Each exercise is three or four cells: the question as
                        the workbook states it, an empty cell pre-seeded with
                        the imports the solution uses, the answer, and for
                        computational ones a cell that recomputes it.
                    </p>
                </div>
                <div class="panel">
                    <h3>Import, modify, re-run</h3>
                    <pre
                        class="small"
                        style="text-align: left; margin: 0.15em 0 0"
                    ><code>from exercises.exercise_09_1 import ANSWER, solve

show(ANSWER)        # the worked prose
result = solve()    # recompute it yourself</code></pre>
                    <p class="small" style="margin-top: 0.35em">
                        <code>solve()</code> takes no arguments and reads the
                        data you already downloaded. Copy its body into your
                        cell, change it, and compare; that is the exercise.
                    </p>
                </div>
            </div>
            <p class="small center muted" style="margin-top: 0.4em">
                The decks are <strong>generated</strong>. Edit the module or the
                chapter, then rebuild in the same container:
                <code
                    >docker compose run --rm jupyter uv run python
                    scripts/make_exercise_notebooks.py</code
                >.
            </p>
            <aside class="notes">
                (~2 min) Walk the shape once, because it is the part people get
                wrong. Everything is named after the workbook: chapter 9
                exercise 1 is exercise_09_1.py, so there is never a question
                about which module answers what. Every module exposes exactly
                two things, ANSWER and solve(), so once they have seen one they
                have seen all 56. The scratch cell is worth pointing at on
                screen: the imports in it are read out of the solution module
                automatically, so it tells them which part of the codebase
                already does the work without handing over the answer. Warn
                them that the notebooks are generated from the modules and the
                chapter text, so an edit to the .ipynb is destroyed by the next
                rebuild and CI's drift check will catch it first. Their PR
                should change a module or a chapter, then include the rebuilt
                notebooks. Everything here runs in the container, including the
                rebuild, so nobody needs Python or uv on their laptop; the
                mounted notebooks folder means the cells they edit are saved
                into their own checkout and go into the PR. Last thing:
                disagreeing with one of our answers is a
                welcome PR. Every number in those answers was run, not guessed,
                but they were run by us, and chapter 9 is the whole argument for
                checking rather than trusting.
            </aside>
        </section>

        <!-- 50 · Hands-on -->
        <section class="dense title-slide center">
            <div class="eyebrow">Hands-on module</div>
            <h1>Run it yourself</h1>
            <p class="subtitle">
                Reproduce and extend Kos et al. 2017 (t-SNE vs UMAP vs EVoC) in
                ~1 minute. Everything runs in Docker; nothing to install.
            </p>
            <pre
                class="small"
                style="text-align: left; max-width: 30em; margin: 0.6em auto"
            ><code>cd day_4_clustering
./run.sh download --all
./run.sh run --fast
docker compose up      # localhost:9999</code></pre>
            <p class="byline">
                <a
                    href="https://github.com/garciadias/iaa-advanced-neural-networks-2026-draft"
                    target="_blank"
                    rel="noopener"
                    >github.com/garciadias/iaa-advanced-neural-networks-2026-draft</a
                >
            </p>
            <p class="small muted" style="margin-top: 0.35em">
                Full walkthrough &middot;
                <code>docs/student_activities.md</code> &nbsp;&middot;&nbsp;
                cluster assignment &middot;
                <code>docs/cluster_assignment.md</code> &nbsp;&middot;&nbsp;
                container details &middot; <code>docs/docker.md</code>
            </p>
            <div
                class="cols"
                style="
                    --n: 3;
                    gap: 0.7em;
                    max-width: 22em;
                    margin: 0.8em auto 0;
                    align-items: stretch;
                "
            >
                <div>
                    <p
                        class="figure-placeholder__desc"
                        style="margin-bottom: 0.35em"
                    >
                        Workshop repo
                    </p>
                    <div
                        class="figure"
                        style="
                            display: block;
                            aspect-ratio: 1 / 1;
                            padding: 0.4em;
                        "
                    >
                        <img
                            :src="asset('presentation.png')"
                            alt="QR code linking to the workshop repository and these slides"
                            title="Workshop repository"
                            style="
                                display: block;
                                width: 100%;
                                height: 100%;
                                object-fit: contain;
                            "
                        />
                    </div>
                </div>
                <div>
                    <p
                        class="figure-placeholder__desc"
                        style="margin-bottom: 0.35em"
                    >
                        Flip a flag
                    </p>
                    <div
                        class="panel"
                        style="
                            aspect-ratio: 1 / 1;
                            min-height: 0;
                            padding: 0.55em;
                            display: flex;
                            align-items: center;
                            justify-content: center;
                        "
                    >
                        <code
                            style="
                                font-size: 0.5em;
                                text-align: left;
                                line-height: 1.5;
                            "
                            >--fast<br />--full<br />--cluster "M 67"<br />--region
                            30</code
                        >
                    </div>
                </div>
                <div>
                    <p
                        class="figure-placeholder__desc"
                        style="margin-bottom: 0.35em"
                    >
                        Three methods
                    </p>
                    <div
                        class="panel flip"
                        style="
                            aspect-ratio: 1 / 1;
                            min-height: 0;
                            padding: 0.55em;
                            display: flex;
                            align-items: center;
                            justify-content: center;
                        "
                    >
                        <span style="font-size: 0.75em; text-align: center"
                            >t-SNE · UMAP · EVoC</span
                        >
                    </div>
                </div>
            </div>
            <aside class="notes">
                (~2 min) Close on the workshop. Everything runs in Docker, so
                the only prerequisites are the ones from section B of the
                school install guide; there is no pip, no conda and no Python
                version to argue with. Three steps: pull the DR19 Astra ASPCAP
                file (1.17 GB, so do it on the hotel wifi tonight, not now),
                one run, one notebook server. run.sh is only a wrapper that
                types the mount flags for them; docker compose up is the same
                image serving JupyterLab on 9999. The fast run caps the field
                at 25 000 stars and finishes in about a minute; `--full` drops
                the cap and takes ten to twenty. Invite them to add
                `--cluster "M 67" --region 30` and watch the precision column
                jump; that is lesson one from the previous slide, reproduced on
                their own laptop in sixty seconds.
            </aside>
        </section>

        <!-- 50b · Where this sits in the literature -->
        <section class="denser">
            <div class="eyebrow">
                Positioning &middot; the honest comparison
            </div>
            <h2>Two papers you should ask me about</h2>
            <div class="cols" style="--n: 2; margin-top: 0.4em">
                <div class="panel">
                    <h3>
                        Hogg et al. 2016, &ldquo;Chemical tagging
                        <em>can</em> work&rdquo;
                    </h3>
                    <p class="small">
                        K-means on 15 <em>Cannon</em> abundances for ~10<sup
                            >5</sup
                        >
                        APOGEE stars, no positional information, and the
                        abundance-space overdensities <em>are</em> phase-space
                        clusters. The strongest published case that abundances
                        alone suffice.
                    </p>
                    <p class="small muted" style="margin-top: 0.3em">
                        Our claim is narrower than &ldquo;abundances
                        fail&rdquo;: on <em>our</em> 16 ASPCAP elements, on
                        <em>these</em> clusters, the masked latent separates
                        better. Hogg's precision (~0.04 dex, Cannon) is not our
                        precision.
                    </p>
                </div>
                <div class="panel flip">
                    <h3>Spina et al. 2025, deep chemical tagging</h3>
                    <p class="small">
                        A graph-attention autoencoder over chemistry + orbits +
                        age, ~47 000 APOGEE thin-disk stars &rarr; 282 groups,
                        recovering 5 of 6 open clusters plus known moving
                        groups.
                        <span class="muted">A&amp;A 702, A267</span>
                    </p>
                    <p class="small muted" style="margin-top: 0.3em">
                        The closest competitor, and
                        <strong>we have not benchmarked against it</strong>.
                        They inject kinematics and age; we deliberately do not.
                        Different question, overlapping claim, the honest next
                        experiment.
                    </p>
                    <p class="small" style="margin-top: 0.3em">
                        <strong>And they tested our idea.</strong> &sect;5.3:
                        clustering their 4-D autoencoder
                        <em>latent</em> recovers only
                        <strong>3 of 6</strong> open clusters, against
                        <strong>5 of 6</strong> in the reconstructed 10-D output
                        and <strong>0 of 6</strong>
                        in the raw abundances. Their reason: a latent
                        &ldquo;loses fine chemical details&rdquo; and &ldquo;is
                        not regularized to be continuous or well-structured for
                        clustering&rdquo;.
                    </p>
                </div>
            </div>
            <p class="small muted center" style="margin-top: 0.45em">
                What is genuinely new here: the features come from the
                <strong>raw spectrum with no labels at all</strong>, not Cannon
                abundances, not ASPCAP abundances, not orbits or ages.
            </p>
            <aside class="notes">
                (~1.5 min) Put this in before someone in the audience does. Hogg
                2016 is titled "Chemical Tagging Can Work" and it is the
                counter-result to the framing of this talk, so state it yourself
                and scope the claim: we are not saying abundances cannot tag, we
                are saying that on our ASPCAP precision and our cluster set, the
                self-supervised latent tags better. Spina 2025 is the direct
                competitor, published last year, doing deep chemical tagging
                with graph attention networks; we have not benchmarked against
                them, and the right answer to "why not" is "not yet, and it's
                the obvious next run", not hand-waving. What survives both
                comparisons is the labels-free part: our features never saw an
                element ratio, a velocity or an age. One more thing you must be
                ready for, because it is the sharpest question in the room:
                Spina's section 5.3 ran the experiment we are advocating
                (cluster the autoencoder's latent) and it did worse than their
                reconstructed output, 3 of 6 clusters against 5 of 6. Do not
                hide it; it is on the slide. The answer has three parts. One,
                their latent is 4-D compressing 10 abundances; ours is 256-D
                compressing 8575 pixels, so "the bottleneck discards fine
                detail" bites very differently. Two, their input is already
                ASPCAP's lossy summary, so their autoencoder compresses a
                compression; ours starts from the raw spectrum (so does our PCA
                baseline, which is why PCA, not the abundances, is the
                comparison that matters). Three, and most honestly: their
                finding is a real warning that reconstruction-trained latents
                are not optimised for clustering, and it is exactly why our
                head-to-head against PCA is close. If someone pushes, concede
                that a like-for-like test, cluster our decoder output as well as
                our latent, is a run we have not done and should.
            </aside>
        </section>

        <!-- 51 · References -->
        <section class="denser">
            <div class="eyebrow">References</div>
            <h2>Papers &amp; code behind this talk</h2>
            <div
                class="cols"
                style="--n: 3; margin-top: 0.3em; align-items: start"
            >
                <div class="panel">
                    <h3>The science &amp; this work</h3>
                    <ul
                        style="
                            list-style: none;
                            margin: 0;
                            font-size: 0.5em;
                            line-height: 1.1;
                        "
                    >
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://arxiv.org/abs/1709.00794"
                                target="_blank"
                                rel="noopener"
                                >Kos et al. 2017</a
                            >, GALAH: chemical tagging of star clusters &amp;
                            new members in the Pleiades.
                            <span class="muted">arXiv:1709.00794</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://ui.adsabs.harvard.edu/abs/2002ARA%26A..40..487F/abstract"
                                target="_blank"
                                rel="noopener"
                                >Freeman &amp; Bland-Hawthorn 2002</a
                            >, The New Galaxy: signatures of its formation.
                            <span class="muted">ARA&amp;A 40, 487</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.3847/1538-4357/833/2/262"
                                target="_blank"
                                rel="noopener"
                                >Hogg et al. 2016</a
                            >, Chemical tagging <em>can</em> work: phase-space
                            structures found by abundance similarity alone.
                            <span class="muted">ApJ 833, 262</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1051/0004-6361/201732134"
                                target="_blank"
                                rel="noopener"
                                >Garcia-Dias et al. 2018</a
                            >, Machine learning in APOGEE: unsupervised spectral
                            classification with K-means.
                            <span class="muted">A&amp;A 612, A98</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1051/0004-6361/201935223"
                                target="_blank"
                                rel="noopener"
                                >Garcia-Dias et al. 2019</a
                            >, Machine learning in APOGEE: stellar populations.
                            <span class="muted">A&amp;A 629, A34</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1016/B978-0-12-815739-8.00013-4"
                                target="_blank"
                                rel="noopener"
                                >Garcia-Dias et al. 2020</a
                            >, Clustering analysis (book chapter).
                            <span class="muted">Elsevier</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1051/0004-6361/202141779"
                                target="_blank"
                                rel="noopener"
                                >Casamiquela et al. 2021</a
                            >, The (im)possibility of strong chemical tagging.
                            <span class="muted">A&amp;A 654, A151</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1093/mnras/staa2743"
                                target="_blank"
                                rel="noopener"
                                >Kreckel et al. 2020</a
                            >, Measuring the mixing scale of the ISM within
                            nearby spiral galaxies.
                            <span class="muted">MNRAS 499, 193</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1051/0004-6361/202555794"
                                target="_blank"
                                rel="noopener"
                                >Spina et al. 2025</a
                            >, Deep chemical tagging: open clusters &amp; moving
                            groups with graph attention networks.
                            <span class="muted">A&amp;A 702, A267</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://arxiv.org/abs/2512.16558"
                                target="_blank"
                                rel="noopener"
                                >Bot, McInnes &amp; Aerts 2025</a
                            >, Persistent multiscale density-based clustering
                            (PLSCAN).
                            <span class="muted">arXiv:2512.16558</span>
                        </li>
                    </ul>
                </div>
                <div class="panel">
                    <h3>Clustering &amp; density algorithms</h3>
                    <ul
                        style="
                            list-style: none;
                            margin: 0;
                            font-size: 0.5em;
                            line-height: 1.1;
                        "
                    >
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1109/TIT.1982.1056489"
                                target="_blank"
                                rel="noopener"
                                >Lloyd 1982</a
                            >, Least squares quantization in PCM (K-means).
                            <span class="muted">IEEE TIT</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            MacQueen 1967, Some methods for classification &amp;
                            analysis of multivariate observations (K-means).
                            <span class="muted">Berkeley Symp.</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf"
                                target="_blank"
                                rel="noopener"
                                >Arthur &amp; Vassilvitskii 2007</a
                            >
                            k-means++: the advantages of careful seeding.
                            <span class="muted">SODA</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://www.aaai.org/Papers/KDD/1996/KDD96-037.pdf"
                                target="_blank"
                                rel="noopener"
                                >Ester et al. 1996</a
                            >, A density-based algorithm for discovering
                            clusters (DBSCAN). <span class="muted">KDD</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            Ankerst et al. 1999, OPTICS: ordering points to
                            identify the clustering structure.
                            <span class="muted">SIGMOD</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1007/978-3-642-37456-2_14"
                                target="_blank"
                                rel="noopener"
                                >Campello et al. 2013</a
                            >, Density-based clustering via hierarchical density
                            estimates. <span class="muted">PAKDD</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1145/2733381"
                                target="_blank"
                                rel="noopener"
                                >Campello et al. 2015</a
                            >, Hierarchical density estimates (HDBSCAN*).
                            <span class="muted">ACM TKDD</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://arxiv.org/abs/1705.07321"
                                target="_blank"
                                rel="noopener"
                                >McInnes &amp; Healy 2017</a
                            >, Accelerated hierarchical density clustering (the
                            <code>hdbscan</code> library).
                            <span class="muted">arXiv:1705.07321</span>
                        </li>
                    </ul>
                </div>
                <div class="panel">
                    <h3>Embeddings &amp; validation</h3>
                    <ul
                        style="
                            list-style: none;
                            margin: 0;
                            font-size: 0.5em;
                            line-height: 1.1;
                        "
                    >
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf"
                                target="_blank"
                                rel="noopener"
                                >van der Maaten &amp; Hinton 2008</a
                            >, Visualizing data using t-SNE.
                            <span class="muted">JMLR 9</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            van der Maaten 2014, Accelerating t-SNE using
                            tree-based algorithms (Barnes-Hut).
                            <span class="muted">JMLR 15</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://arxiv.org/abs/1802.03426"
                                target="_blank"
                                rel="noopener"
                                >McInnes, Healy &amp; Melville 2018</a
                            >, UMAP: uniform manifold approximation &amp;
                            projection.
                            <span class="muted">arXiv:1802.03426</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1016/0377-0427(87)90125-7"
                                target="_blank"
                                rel="noopener"
                                >Rousseeuw 1987</a
                            >, Silhouettes: a graphical aid.
                            <span class="muted">J. Comput. Appl. Math.</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1214/aos/1176346577"
                                target="_blank"
                                rel="noopener"
                                >Hartigan &amp; Hartigan 1985</a
                            >, The dip test of unimodality.
                            <span class="muted">Ann. Stat.</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://aclanthology.org/D07-1043.pdf"
                                target="_blank"
                                rel="noopener"
                                >Rosenberg &amp; Hirschberg 2007</a
                            >, V-measure: external cluster evaluation.
                            <span class="muted">EMNLP</span>
                        </li>
                        <li style="margin-bottom: 0.12em">
                            <a
                                href="https://doi.org/10.1016/j.neuroimage.2009.06.014"
                                target="_blank"
                                rel="noopener"
                                >Nanetti et al. 2009</a
                            >, Repeated K-means cortical parcellation.
                            <span class="muted">NeuroImage</span>
                        </li>
                    </ul>
                </div>
            </div>
            <p class="small" style="margin-top: 0.45em">
                Code:
                <a
                    href="https://github.com/TutteInstitute/evoc"
                    target="_blank"
                    rel="noopener"
                    >EVoC</a
                >
                ·
                <a
                    href="https://github.com/jelmerbot/fast_plscan"
                    target="_blank"
                    rel="noopener"
                    >PLSCAN</a
                >
                ·
                <a
                    href="https://github.com/garciadias/iaa-advanced-neural-networks-2026-draft"
                    target="_blank"
                    rel="noopener"
                    >workshop repo</a
                >
            </p>
            <aside class="notes">
                (~1 min) Appendix, do not present line-by-line. Leave it up
                while taking questions, or point people here for the arXiv/DOI
                links. Every claim in the talk traces to one of these. Two
                entries carry no link on purpose: MacQueen 1967 (proceedings, no
                stable DOI) and van der Maaten 2014, where the JMLR volume needs
                checking before a URL goes on the slide.
            </aside>
        </section>
    </RevealDeck>
</template>

<style>
/*
  Deck-wide backdrop: public/static/img/astro_background.png — a spectrum wall
  dissolving into a star field.

  The image is dark exactly where the text lives (42% of its pixels below
  luminance 0.05; the right third is a black star field, the lower left
  blue-violet brick), and this deck's ink is dark by design, so the picture
  cannot sit behind text unmixed: on a flat veil a 13px caption over the dark
  third reads 1.5:1 at opacity 0.55 and needs ~0.95 to reach 4.5:1 — at which
  point there is no picture left. So the image is a *mat* and the slide itself
  is the readable surface:

  - .reveal carries the image (veiled to 0.5) and shows it in the frame around
    the canvas — 20-35px in present mode, the page margin on the site.
  - .slides > section is the slide surface, near-opaque, so everything the deck
    draws sits on it: bare eyebrows, headings, muted captions, both tables, the
    About-me timeline's SVG labels, the translucent chips. 0.94 keeps a 6% tint
    of the picture inside the canvas and holds the ink at full contrast.

  Measured over all 86 slides before this: the deck's own floor was 4.32:1 (80
  of 903 text elements, every one a muted-grey caption on --surface). Over this
  surface, measured across 4,770 ink samples, the floor is 4.56:1 and the muted
  copy sits at 5.7–6.1:1 — nothing in the deck is below AA any more.

  Dials: the 0.5 on the image (frame strength) and the 0.94 on the surface.
  Lowering the surface past ~0.9 starts letting the star field into text boxes.
*/
.reveal.deck-theme.astro-theme,
html:not(.dark) .reveal.deck-theme.astro-theme,
html.dark .reveal.deck-theme.astro-theme {
    /* Both mode selectors are spelled out on purpose: the theme's own rules are
     scoped `html:not(.dark) .reveal.deck-theme` / `html.dark .reveal.deck-theme`,
     which carry one more simple selector than `.reveal.deck-theme.astro-theme`
     and therefore beat every custom property the astro palette sets — which is
     why this deck has been rendering with the base light palette (and its
     #64748b ink) rather than the astro one it asks for. Matching their
     specificity here is what makes the backdrop apply. */
    background-image:
        linear-gradient(rgba(251, 251, 253, 0.5), rgba(251, 251, 253, 0.5)),
        url("/static/img/astro_background.png");
    background-size: cover, cover;
    background-position: center, center;
    background-repeat: no-repeat, no-repeat;
}

/* The slide surface: opaque enough that no text in the deck depends on what is
   behind it, slightly transparent so the picture still tints the canvas.

   On .slides rather than on the sections: a section is only as tall as its own
   content (336px for the title slide), so a surface there paints a band and
   leaves the image above and below it. .slides *is* the canvas.

   Rounded to 12px, matching the deck's own card radius (.figure and
   .figure-placeholder are 12px, .panel and .slide-body 14px). The radius is not
   just taste: on 54 slides a card is flush with the canvas edge — 21 figures, 26
   panels, 5 slide-body blocks — so the cut has to be no larger than the smallest
   of them. At 12px the flush figures' corners sit concentric with the surface's;
   at 18px the cut went deeper than their own corners and left a 4px sliver of
   white card sitting against the image. It is paint-only — no overflow
   clipping — so no content moves or gets cut, and the one full-bleed slide (the
   video) takes the same radius explicitly below.

   The surface is drawn 3px *outside* the canvas on every side so the content
   never sits on its edge — grown outward rather than insetting the content,
   which is what the flush cards would otherwise feel. That means one fill on a
   pseudo-element rather than a background plus a box-shadow of the same colour:
   two separately antialiased edges at the canvas boundary leave a 1px hairline
   of image showing through (measured: rgb(249,227,223) against an interior of
   rgb(251,245,246)), where a single fill has only its outer edge to antialias.
   Radius 15 = 12 + 3 about the same centre, so the arc stays concentric with the
   flush cards' corners. */
.reveal.deck-theme.astro-theme .slides,
html:not(.dark) .reveal.deck-theme.astro-theme .slides,
html.dark .reveal.deck-theme.astro-theme .slides {
    background: none;
}

.reveal.deck-theme.astro-theme .slides::before,
html:not(.dark) .reveal.deck-theme.astro-theme .slides::before,
html.dark .reveal.deck-theme.astro-theme .slides::before {
    content: "";
    position: absolute;
    inset: -3px;
    background: rgba(251, 251, 253, 0.94);
    border-radius: 15px;
}

/* The video slide is drawn edge to edge by YouTube's own iframe, so it would
   square off the corners every other slide rounds: give it the same radius and
   clip, since a full-bleed element is the one case where the surface's paint-only
   rounding needs overflow to follow it. No mode prefixes needed here — nothing
   else styles this class. */
.reveal.deck-theme.astro-theme section.video-slide {
    border-radius: 12px;
    overflow: hidden;
}

/* The branding chrome (logo + QR chip) sits below the slides on purpose — slide
   content must be able to cover it — but that puts it below the surface above
   as well, where it would disappear. Raised for this deck so the lecture QR
   stays scannable on every slide. */
.reveal.deck-theme.astro-theme .deck-brand-bottom,
html:not(.dark) .reveal.deck-theme.astro-theme .deck-brand-bottom,
html.dark .reveal.deck-theme.astro-theme .deck-brand-bottom {
    z-index: 2;
}

/* On the slides whose content already occupies the chip's corner the chip steps
   aside (measured in the script block); everywhere else it stays above the
   surface, which is what keeps the QR scannable in the lecture. */
.reveal.deck-theme.astro-theme .deck-brand-bottom.qr-chip-covered,
html:not(.dark)
    .reveal.deck-theme.astro-theme
    .deck-brand-bottom.qr-chip-covered,
html.dark .reveal.deck-theme.astro-theme .deck-brand-bottom.qr-chip-covered {
    display: none;
}
</style>
