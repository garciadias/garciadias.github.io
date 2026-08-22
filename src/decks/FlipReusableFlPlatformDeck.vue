<script setup>
import DeckBrand from '@/components/DeckBrand.vue';
import RevealDeck from '@/components/RevealDeck.vue';

// 25-minute talk for "From Zero to Hero: Federated AI in Healthcare Systems"
// (St John's College, Cambridge). Audience: clinicians, researchers and people
// new to federated learning — so the argument is deliberately platform-first
// rather than algorithm-first.
//
// Narrative, in four movements (slide numbers in the comments below match the
// order of the sections in this file, 0-19):
//   0-3   the problem: the data exists, cannot move, and its diversity is the
//         point → Flower and NVIDIA FLARE solved the federated core, so what is
//         left unsolved is the layer between a hospital and the framework.
//   4-11  the platform: physical architecture, the services it decomposes into,
//         the network ask, where the guardrails execute, deployment modes,
//         governance, and the job-type contract that makes it reusable.
//   12-17 the demo: a real UK ⇄ Thailand deployment (the organisers asked for
//         this specifically, so it sits at the centre of the talk, ~7 min) —
//         the four researcher steps, per-site approval, a live run, then two
//         recorded walkthroughs on unrelated workloads.
//   18-19 the team, and the close.
//
// deckAsset()  → this deck's own folder: bridge artwork, the intro globe, the
//                UI recordings, the two demo-video QR codes.
// asset()      → this deck reuses the flip-maturity-pitch-2026 image folder
//                rather than duplicating ~4 MB of screenshots.
// sharedAsset() → public/presentations/ (team headshots).
// pydataAsset() → the QR codes built for the PyData lightning talk.
const deckAsset = (name) => `${import.meta.env.BASE_URL}presentations/flip-reusable-fl-platform/${name}`
const asset = (name) => `${import.meta.env.BASE_URL}presentations/flip-maturity-pitch-2026/${name}`
const sharedAsset = (name) => `${import.meta.env.BASE_URL}presentations/${name}`
const pydataAsset = (name) => `${import.meta.env.BASE_URL}presentations/pydata-london-flip-lightning-2026/${name}`
</script>

<template>
  <!-- --hero-poster rides on the .reveal root so it reaches reveal's own
       .backgrounds layer, which lives outside this component's markup. It paints
       the intro globe's still frame behind the video: reveal has no poster
       option for background video, and a browser that blocks autoplay would
       otherwise show a blank white slide. -->
  <RevealDeck
    :options="{ center: true }"
    :style="{ '--hero-poster': `url(${deckAsset('intro-globe-poster.jpg')})` }"
  >
    <template #chrome>
      <DeckBrand :qr="deckAsset('presentation.png')" />
    </template>

    <!-- 0 · Title. The animated globe is lifted from the FLIP_2026_08 deck:
         reveal's own data-background-video layer handles the full-bleed sizing
         and only plays while this slide is active. The poster covers the first
         paint and any browser that refuses to autoplay. -->
    <section
      class="title-slide center video-hero"
      :data-background-video="`${deckAsset('intro-globe.mp4')},${deckAsset('intro-globe.webm')}`"
      data-background-video-loop
      data-background-video-muted
      data-background-size="cover"
      data-background-color="#ffffff"
    >
      <div class="eyebrow">From Zero to Hero: Federated AI in Healthcare Systems · St John's College, Cambridge</div>
      <h1>FLIP: a reusable, general-purpose open-source FL platform</h1>
      <div class="slide-body">
        <p class="byline">
          Rafael Garcia-Dias
        </p>
        <p class="venue-note">github.com/londonaicentre/FLIP · Apache 2.0</p>
      </div>
      <aside class="notes">
        (~45s) Open with the thesis, not the CV. Most FL talks are about algorithms; this one is about
        the boring 90% that decides whether an FL collaboration happens at all. FLIP is our answer:
        one platform, installed once per site, reused across projects — open source, Apache 2.0.
        Everything in this talk is backed by a running deployment and dated screenshots, not a roadmap.
        Say up front there is a live demo in the middle, so the room knows to stay awake for it.
        The globe loops for ~7s — let it run while you open, it earns its keep as the room settles.
      </aside>
    </section>

    <!-- 1 · Roadmap. Four beats, so the room knows the demo is coming and can
         hold its questions until the end. -->
    <section>
      <div class="eyebrow">What you'll get in the next 25 minutes</div>
      <h2>FLIP — Federated Learning Interoperability Platform</h2>
      <div class="slide-body medium">
        <ol class="contribs">
          <li><strong>The gap</strong> Why are the FL frameworks ready, but not yet mainstream in hospitals?</li>
          <li><strong>The platform</strong> What FLIP is, what it's built on, and what it deliberately doesn't do?</li>
          <li><strong>The demo</strong> A cross-continental deployment, cohort query, governance, training run, results.</li>
          <li><strong>What's next</strong> How to join the community and contribute?</li>
        </ol>
      </div>
      <aside class="notes">
        (~45s) Set expectations explicitly, especially the last line: this audience will hear plenty of
        algorithm talks over two days. Position this one as the deployment counterweight. Signpost that
        the demo is roughly the middle third of the talk.
      </aside>
    </section>

    <!-- 2 · The problem, 1 of 2. Open on scale, not scarcity, then say why the
         data is unreachable. The NHS imaging chart is reused from the PyData
         lightning talk (hence pydataAsset) rather than copied into this deck's
         folder. -->
    <section>
      <div class="eyebrow">The problem · 1 of 2</div>
      <h2>Unlocking the potential of medical imaging data</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1.3fr 0.95fr; margin-top: 0.2em; align-items: start">
          <div>
            <p class="small">
              Every major NHS trust has a PACS archive going back decades with annotated X-rays, CTs and MRIs.
              And yet almost none of it has ever trained a model at national scale.
            </p>
            <ul class="crosslist small">
              <li><strong>Data is fragmented</strong> It sits locked inside individual hospitals and health systems.</li>
              <li><strong>It can't just move</strong> Privacy law, governance processes and cross-border rules mean patient data cannot simply be pooled into a data lake.</li>
              <li><strong>The diversity is the point</strong> A model trained on one country's scanners, protocols and population is exactly the model that fails somewhere else.</li>
            </ul>
          </div>
          <div class="figure" style="aspect-ratio: 4 / 3">
            <img
              :src="pydataAsset('nhs_medical_images.png')"
              alt="Line chart of annual NHS medical imaging procedures in England, rising from 39.94M in 2014/15 to 47.16M in 2023/24, with a dip to 34.92M during the 2020/21 COVID lockdown"
            />
          </div>
        </div>
      </div>
      <aside class="notes">
        (~45s) Lead with scale, not scarcity: this is not a small-data problem, it is a huge,
        already-labelled dataset nobody can train on yet — the chart is England alone, ~47M imaging
        procedures a year and climbing. Then turn to why it is unreachable: healthcare data is siloed
        by design, for good reasons. The third bullet is the one worth landing for a clinical
        audience — this isn't only about volume, it's about generalisation across scanners,
        protocols and populations.
      </aside>
    </section>

    <!-- 3 · The problem, 2 of 2 — the thesis slide, and the load-bearing one.
         Framed as "the frameworks solved their half, and solved it well; this is
         the half nobody owns yet". Flower Labs sponsors this event and has
         committed code to FLIP, so the argument builds on the FL stacks rather
         than positioning FLIP against them. The pills here are the agenda for
         slides 6-11, and the two piers set up the two backends named on 6. -->
    <section>
      <div class="eyebrow">The problem · 2 of 2 · the thesis of this talk</div>
      <h2>On the shoulders of giants</h2>
      <div class="slide-body">
        <!-- align-items: start keeps the two columns independent: the text sets its own
             height and stays anchored under the heading instead of being pushed down to
             the vertical centre of whatever the (taller) figure column happens to be.
             flip_bridge.png is 992×1074 and nearly square, so it is wide enough to force
             the grid columns open on its own — min-width: 0 lets the cell shrink, and the
             image is capped in both axes so the ratio below is what actually decides the
             split. -->
        <div class="fig-split" style="--cols: 1.65fr 1fr; margin-top: 0; align-items: start">
          <div>
            <div class="panel flip">
              <h3>The federated core is ready</h3>
              <p class="small">
                <strong>Flower</strong> and <strong>NVIDIA FLARE</strong> made it something you can
                build and depend on. Aggregation strategies, client/server orchestration, secure transport,
                real communities behind both.
              </p>
            </div>
            <p class="small" style="margin-top: 0.5em">
              The bridge <em>between a hospital and the FL network</em>
            </p>
            <div style="margin-top: 0.35em">
              <span class="pill">cohort definition</span><span class="pill">data harmonisation</span
              ><span class="pill">imaging retrieval</span><span class="pill">per-site approval</span
              ><span class="pill">governance</span><span class="pill">job scheduling</span
              ><span class="pill">monitoring</span><span class="pill">audit trail</span
              ><span class="pill">networking</span><span class="pill">provenance</span
              ><span class="pill">permission controls</span>
            </div>
          </div>
          <div class="figure" style="justify-self: center; min-width: 0; max-width: 100%">
            <img
              :src="deckAsset('flip_bridge.png')"
              alt="Illustration: two hospitals on opposite sides of a river valley, each with its own PACS, EHR and governance office. A stone bridge labelled FLIP — Federated Learning Interoperability Platform — spans the gap between them, carried on two piers marked Flower and NVIDIA FLARE"
              style="max-width: 100%; max-height: 470px"
            />
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min) The load-bearing slide of the talk — and the one where tone matters most, with
        Flower in the room as a sponsor. Walk the drawing: the two piers are Flower and NVIDIA
        FLARE, and everything above them only stands because they are there. We depend on both,
        and Flower Labs engineers have committed code to FLIP. The gap this talk is about is the
        deck of the bridge — from a hospital's own PACS, EHR and governance office out to the FL
        network. It isn't anyone's research contribution, so it gets rebuilt by hand every time.
        Admit we were guilty of it too: the first version of FLIP was project-shaped. The pivot
        was deciding the platform outlives the project.
      </aside>
    </section>

    <!-- 4 · Architecture, 1 of 2 — the symmetric picture plus the four
         properties a site cares about. Deliberately shares its heading with the
         next slide: same claim, then the same claim in detail. -->
    <section>
      <div class="eyebrow">The platform · how it physically works</div>
      <h2>One Central Hub, one Trust Node per site</h2>
      <div class="slide-body">
        <div style="display: grid; grid-template-rows: auto auto; gap: 0.5em; align-items: center">
          <div class="figure center">
            <img
              :src="asset('flip-architecture-symmetric.png')"
              alt="Symmetric FLIP architecture: two mirrored trust nodes, each taking PACS, EHR and other clinical systems into a secure enclave holding a structured OMOP database and XNAT, feeding an AI training client, connected through a FLIP node to the platform"
              style="max-height: 430px; max-width: 100%"
            />
          </div>
          <ul class="dotlist small" style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 0.8em; margin: 0; padding: 0">
            <li><strong>Central Hub</strong><br>Coordinates jobs, holds no patient data. Cloud-hosted or on-prem.</li>
            <li><strong>Trust Node</strong><br>Lives inside the hospital's own secure enclave: an OMOP database, an XNAT imaging archive fed from PACS, and the FL client on local GPU compute.</li>
            <li><strong>Outbound-only</strong><br>The node calls out over HTTPS. No inbound firewall rule is ever requested.</li>
            <li><strong>Only weights, aggregate statistics and log messages cross the wire</strong><br>Images and row-level records never leave the building.</li>
          </ul>
        </div>
      </div>
      <aside class="notes">
        (~1 min) Walk the diagram right to left for a clinical audience: PACS and the EHR are systems
        they already know; the enclave is the bit their IT department already runs; FLIP is the thin
        layer that makes those two things addressable by a researcher who will never log in. The
        outbound-only point is the single most effective thing to say to a hospital network team —
        approving one outbound HTTPS connection is a categorically easier conversation than an
        inbound rule.
      </aside>
    </section>

    <!-- 5 · Architecture, 2 of 2 — the full component diagram, given the whole
         slide. Same eyebrow, heading and notes as slide 4 on purpose; if only
         one of the two survives a timing cut, either works alone. -->
    <section>
      <div class="eyebrow">The platform · how it physically works</div>
      <h2>One Central Hub, one Trust Node per site</h2>
      <div class="slide-body">
        <div class="figure center" style="width: 100%">
          <img
            :src="deckAsset('flip_architecture-flip_architecture.png')"
            alt="Detailed FLIP architecture: a cloud-hosted FLIP Hub holding the FLIP UI, Hub API, FL API and FL Server; on the right, a hospital's EHR and PACS feeding a structured OMOP database and an XNAT imaging archive inside a secure enclave, each read by its own Data API and Imaging API, with the Trust API orchestrating both. Only two links cross the FLIP node boundary: Hub API to Trust API, and FL Server to FL Client"
            style="max-height: 640px; max-width: 100%"
          />
        </div>
      </div>
      <aside class="notes">
        (~1 min) Walk the diagram right to left for a clinical audience: PACS and the EHR are systems
        they already know; the enclave is the bit their IT department already runs; FLIP is the thin
        layer that makes those two things addressable by a researcher who will never log in. The
        outbound-only point is the single most effective thing to say to a hospital network team —
        approving one outbound HTTPS connection is a categorically easier conversation than an
        inbound rule.
      </aside>
    </section>

    <!-- 6 · Microservice inventory. Every name here is the real Compose service
         key from the FLIP monorepo, so anyone who clones the repo reads the same
         words. The two "only" claims are the point of the slide: one service
         reaches the OMOP database, one reaches the imaging archive. -->
    <section>
      <div class="eyebrow">The platform · what the boxes on that diagram actually are</div>
      <h2>Small services, one job each</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1fr 1.15fr; margin-top: 0; align-items: start">
          <div class="panel">
            <h3>Central Hub</h3>
            <ul class="dotlist small" style="margin-top: 0.2em">
              <li><code>flip-api</code> projects, cohorts, approvals, scheduling, user management.</li>
              <li><code>flip-ui</code> Vue 3, the web interface.</li>
              <li><code>fl-api-net-N</code> + <code>fl-server-net-N</code> one pair <em>per FL network</em>.</li>
            </ul>
          </div>
          <div class="panel flip">
            <h3>Trust Node</h3>
            <ul class="dotlist small" style="margin-top: 0.2em">
              <li><code>trust-api</code> polls jobs from <code>flip-api</code>, drives the <code>data-access-api</code> and <code>imaging-api</code>, holds the site's keys.</li>
              <li><code>data-access-api</code> the <em>only</em> service that speaks SQL to the anonymised OMOP database.</li>
              <li><code>imaging-api</code> the <em>only</em> service that speaks to XNAT and Orthanc/PACS.</li>
              <li><code>fl-client-net-N</code> on local GPUs, plus a Loki/Grafana stack so logs stay on site.</li>
            </ul>
          </div>
        </div>
        <div class="cols compact" style="--n: 4; margin-top: 0.4em">
          <div class="panel center small">
            <h3>One monorepo</h3>
            <p class="small muted">One command deploys the whole stack locally: <code>make&nbsp;up</code>.</p>
          </div>
          <div class="panel center small">
            <h3>Simple integration</h3>
            <p class="small muted">The <code>flip</code> library lets researchers run local experiments.</p>
          </div>
          <div class="panel center small">
            <h3>Two backends</h3>
            <p class="small muted">NVFLARE and Flower.</p>
          </div>
          <div class="panel center small">
            <h3>Apps</h3>
            <p class="small muted">Customisable modules that extend the platform's functionality.</p>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~45s) The point of this slide is that the architecture diagram is not aspirational — these
        are the container names. Land the two "only" claims: exactly one service can reach the OMOP
        database and exactly one can reach the imaging archive, which is what makes the trust-side
        attack surface reviewable by a hospital security team in an afternoon. The per-net pair is
        the detail people miss: FL servers are not a single shared bottleneck, they are provisioned
        per network. If time is short, read the four panels and move on.
      </aside>
    </section>

    <!-- 7 · The outbound-only mechanism, with the constants taken from the code
         (POLL_INTERVAL_SECONDS = 5, HEARTBEAT_TIMEOUT_SECONDS = 30, six TaskType
         values). This is the slide that answers "what did you open on our
         network?" with a number. -->
    <section>
      <div class="eyebrow">The platform · what we asked the network team for</div>
      <h2>Simple and secure design</h2>
      <div class="slide-body">
        <div class="kpis" style="--n: 4">
          <div class="panel center">
            <div class="stat" style="font-size: 1.5em">0</div>
            <div class="stat-label">inbound ports on the trust host, and no port 22 open either</div>
          </div>
          <div class="panel center">
            <div class="stat" style="font-size: 1.5em">5s</div>
            <!-- Kept as a single code chip: when the endpoint wraps across two
                 lines, its border makes it read as two broken boxes. -->
            <div class="stat-label">poll interval, <code>trust-api</code> asks the hub for <code>/tasks/pending</code></div>
          </div>
          <div class="panel center">
            <div class="stat" style="font-size: 1.5em">30s</div>
            <div class="stat-label">of silence and the hub marks that trust offline in the UI</div>
          </div>
          <div class="panel center">
            <div class="stat" style="font-size: 1.5em">6</div>
            <div class="stat-label">API verbs — the hub's entire request vocabulary</div>
          </div>
        </div>
        <div class="fig-split" style="--cols: 1.25fr 0.8fr; margin-top: 0.5em; align-items: center">
          <ul class="checklist small">
            <li>The hub has <strong>no route into the trust</strong>. Work is a queue the site chooses to read, so a site that pauses simply stops polling.</li>
            <li>The FL clients hold <strong>no Central Hub credentials at all</strong>. Compromising a client doesn't give an attacker access to the hub.</li>
            <li>Registration mints the site's keys and the hub keeps only a <strong>SHA&#8209;256 hash</strong> of the API key; the plaintext lives in the site's own kit file.</li>
          </ul>
          <div class="panel">
            <h3>Trust API verbs</h3>
            <p class="small" style="margin-bottom: 0">
              <code>cohort_query</code> <code>create_imaging</code> <code>get_imaging_status</code>
              <code>delete_imaging</code> <code>reimport_studies</code> <code>update_user_profile</code>
            </p>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min) This is the slide to have ready for anyone from an IT or infosec team, and the
        numbers are the whole argument: nothing to open, a five-second poll outbound over HTTPS, and
        a command vocabulary short enough to print. Say the last bullet plainly — the hub cannot
        reach in, so "unplugging" is not an escalation, it is the site stopping its own poller. The
        six verbs are worth reading out: notice how much of what people fear a platform can do is
        simply not expressible.
      </aside>
    </section>

    <!-- 8 · Where the guardrails physically execute. The asymmetry is deliberate
         and is the interesting engineering claim: the hub's copy of the query
         validator is convenience, the trust's copy is the control. -->
    <section>
      <div class="eyebrow">The platform · governance encoded, not promised</div>
      <h2>Every guardrail runs on the trust's side of the wire</h2>
      <div class="slide-body">
        <div class="cols" style="--n: 3">
          <div class="panel">
            <h3>Parsed, not pattern-matched</h3>
            <p class="small">
              The authoritative validator is the <strong>trust's</strong>: a single statement,
              <code>SELECT</code> only, pinned to the <code>omop</code> schema, a literal row limit,
              re-emitted from the checked syntax tree and run as a <strong>read-only role</strong>.
              No keyword denylist anywhere.
            </p>
          </div>
          <div class="panel">
            <h3>A disclosure floor per site</h3>
            <p class="small">
              Each site sets its own threshold (default <strong>10</strong>). Below it, both
              row-level routes return the <em>same</em> refusal as an empty result, so a researcher
              cannot infer patient-level information by narrowing a query until one person matches.
            </p>
          </div>
          <div class="panel flip">
            <h3>Uploads are quarantined</h3>
            <p class="small">
              Researcher files land in a staging prefix, get structurally scanned, and only then are
              promoted to <code>scanned/</code>. The app bundler reads <strong>only</strong> the
              promoted prefix, so an unscanned file cannot reach a hospital.
            </p>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min) The engineering claim worth landing: the hub and the trust both validate cohort
        queries, and that duplication is deliberate rather than sloppy — the hub's version exists so
        a typo fails in the researcher's hand instead of after a fan-out, while the trust's version
        is the one that is actually load-bearing. The suppression detail in panel two is the one
        governance leads react to: a shared refusal string means a researcher cannot use repeated
        small queries to infer a count. If asked about signature-based AV: structural scanning is in
        place today, signature scanning is tracked separately.
      </aside>
    </section>

    <!-- 9 · Deployment modes — three objections (how, what hardware, what
         network) answered before anyone raises them. -->
    <section>
      <div class="eyebrow">The platform · will it fit our stack?</div>
      <h2>Low-barrier deployment on Trust Nodes</h2>
      <div class="slide-body">
        <div class="cols" style="--n: 3">
          <div class="panel">
            <h3>Deploy how you already operate</h3>
            <p class="small">
              Microservices deployable via <strong>Docker Compose</strong>, a <strong>Kubernetes</strong>
              Helm chart, or <strong>cloud</strong> infrastructure-as-code. A site picks the one its IT team already runs.
            </p>
          </div>
          <div class="panel">
            <h3>Simple hardware requirements</h3>
            <p class="small">
              Minimum to join: a <strong>consumer-grade GPU</strong>, <strong>16&nbsp;GB</strong> RAM,
              <strong>1&nbsp;TB</strong> storage beyond patient data. No datacentre purchase required
              before the first project.
            </p>
          </div>
          <div class="panel flip">
            <h3>Network security teams say yes</h3>
            <p class="small">
              <strong>No ingress connections</strong> to the Trust, ever. All communication is outbound
              and encrypted, compatible with strict hospital network policies.
            </p>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~45s) This slide exists to kill three objections before they're raised. The hardware one is
        the surprise for most people: the entire UK–Thailand validation ran on single workstations
        with one consumer GPU each. Production workloads want more, but the entry ticket is not a
        cluster — it's a machine a research group already has under a desk.
      </aside>
    </section>

    <!-- 10 · Governance as architecture rather than paperwork: the per-site
         veto, the two legal regimes it has already passed, and the guardrails
         that make the veto structural. -->
    <section>
      <div class="eyebrow">The platform · governance is architecture, not paperwork</div>
      <h2>Every site keeps its own veto, on every project</h2>
      <div class="slide-body">
        <div class="cols" style="--n: 3">
          <div class="panel">
            <h3>Sovereign control</h3>
            <p class="small">
              Each institution runs its <strong>own internal governance</strong> for every project.
              Projects are approved or vetoed independently, per site.
            </p>
          </div>
          <div class="panel">
            <h3>Compliance</h3>
            <p class="small">
              Observes the <strong>NHS National Data Opt-Out</strong> in the UK; the Thai node was
              cleared against the <strong>Thailand Personal Data Protection Act</strong> by the hospital
              group's own security contractor.
            </p>
          </div>
          <div class="panel flip">
            <h3>Guardrails</h3>
            <p class="small">
              Researchers only ever see <strong>anonymised</strong> aggregate data at cohort selection.
              Only weights, aggregate statistics and log files leave the site. Every action lands in an
              audit trail.
            </p>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min) The point to land: governance that lives in a signed PDF gets renegotiated every
        project; governance encoded in the platform is enforced by default. The Thailand example is
        the strongest evidence that the model travels — a completely different legal regime, reviewed
        by a third party we had no relationship with, and the answer was still "the node is
        acceptable" because the guarantees are structural.
      </aside>
    </section>

    <!-- 11 · The job-type contract, and the four workload families already
         shipped on it — the evidence that one contract adapts to functionally
         unrelated workflows. -->
    <section>
      <div class="eyebrow">The platform · what it already runs</div>
      <h2>Modular job types, many workflows</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1fr 1.5fr; margin-top: 0; align-items: center">
          <div>
            <p class="small">
              A <strong>job type</strong> fixes the server-side contract: orchestration,
              aggregation and checkpoint staging. The researcher's <strong>app bundle</strong> fills the
              client side. Ship a new app, select a job type, and it runs with
              <strong>no Central Hub or Trust Node code changes</strong>.
            </p>
            <p class="small">
              <em>
              Tutorials:
              <a href="https://github.com/londonaicentre/FLIP/tree/main/fl-tutorials" target="_blank" rel="noopener">londonaicentre/FLIP/fl-tutorials</a>
              </em>
            </p>
          </div>
          <div class="figure center">
            <img
              :src="deckAsset('diagram-components-nets-scheduler.png')"
              alt="Diagram of FLIP's per-run networks: the FLIP Scheduler in the Central Hub brings up one isolated network per training run — net 1 and net 2 — each with its own FastAPI FL API and FL Server in the Central Hub, talking to a GPU-backed FL Client container inside Trust A and Trust B"
              style="width: 100%; height: auto"
            />
            <div class="caption">One isolated network per run, one client per Trust.</div>
          </div>
        </div>
        <div class="cols compact" style="--n: 4; margin-top: 0.4em">
          <div class="panel center">
            <h3>Classification</h3>
            <p class="small muted">CXR classification, Ark+ fine-tuning.</p>
            <span class="pill">FLARE</span><span class="pill">Flower</span>
          </div>
          <div class="panel center">
            <h3>Segmentation</h3>
            <p class="small muted">3D spleen, MONAI bundle.</p>
            <span class="pill">FLARE</span><span class="pill">Flower</span>
          </div>
          <div class="panel center">
            <h3>Evaluation</h3>
            <p class="small muted">Benchmarking, DeLong tests.</p>
            <span class="pill">FLARE</span><span class="pill">Flower</span>
          </div>
          <div class="panel center">
            <h3>Synthesis</h3>
            <p class="small muted">Latent diffusion, federated.</p>
            <span class="pill">FLARE</span>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~45s) Two ideas, one slide. First, the job type is the modularity boundary: the server
        side is fixed per job type, the researcher writes only the client-side app, so a new
        workflow deploys as a bundle + job type with no platform reconfiguration. Second, that
        boundary is not theoretical — the shipped tutorials already span classification,
        segmentation, evaluation and synthesis, each on one or both backends, which is the
        evidence that one contract adapts to functionally unrelated workflows.
      </aside>
    </section>

    <!-- 12 · Demo divider. Deliberate gear change: argument before, product
         on screen after. Carries the demo run-sheet in its notes. -->
    <section class="title-slide center">
      <div class="eyebrow">Part three · roughly a third of this talk</div>
      <h1>Let's look at the real thing</h1>
      <div class="slide-body">
        <p class="subtitle">
          A live cross-continental deployment: <strong>UK ⇄ Thailand</strong>, running today.
        </p>
        <p class="venue-note">Cohort query → governance → training run → results</p>
      </div>
      <aside class="notes">
        (~20s) Deliberate gear change. Everything before this was argument; everything after is the
        product on screen. If running live, this is the moment to switch to the browser — say out loud
        that the screenshots in the deck are the fallback and are from the same system on dated runs.
        DEMO RUN-SHEET (target ~7 min): 1) sign in to the Central Hub; 2) show the two live sites;
        3) run the cohort query and show the per-site aggregate counts; 4) show the approval/DAC
        screen; 5) upload or open an existing app; 6) show the live training feed and loss curves;
        7) show the finished run's provenance record. If the network is unreliable, drive the same
        seven beats from the screenshots and the recorded videos that follow.
      </aside>
    </section>

    <!-- 13 · The deployment being demoed: the estate at both ends, and the task
         it was pointed at, stacked in one figure card. -->
    <section>
      <div class="eyebrow">The demo · the deployment behind it</div>
      <h2>Cross-continental deployment</h2>
      <div class="center">
        <div class="figure center" style="width: 80%">
            <img
              :src="asset('flip-architecture-uk-thailand.png')"
              alt="FLIP architecture: a cloud-hosted Central Hub coordinating a Thailand trust node and a UK trust node, each with EHR/PACS feeding a structured OMOP database and XNAT imaging archive inside a secure enclave"
              style="width: 100%; height: auto"
            />
            <img
              :src="asset('flip-cxr-pathology-classes.png')"
              alt="Five chest X-ray pathology classes used in the study: effusion, consolidation, infiltration, lung nodule or mass, pneumothorax"
              style="width: 100%; height: auto"
            />
            <div class="caption">The shared task: five CXR pathology classes.</div>
        </div>
      </div>
      <aside class="notes">
        (~45s) Why this pairing is the interesting test and not just a nice map: it deliberately
        crosses every axis at once — jurisdiction, sector (public/private), language, network estate,
        and hardware. If the platform abstraction leaked anywhere, it would leak here. Note for
        honesty: this particular study ran on synthetic imaging cohorts. There is no limitations slide
        any more, so say it here rather than being caught by it in Q&amp;A.
      </aside>
    </section>

    <!-- 14 · The researcher's whole job in four steps, with the recorded
         end-to-end GIF covering steps 1-2 beside it. -->
    <section>
      <div class="eyebrow">The demo · the researcher's whole job</div>
      <h2>Four steps, all from a laptop</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1.05fr 1fr; align-items: center">
          <ol class="contribs small">
            <li>
              <strong>Define the cohort.</strong> One SQL query against each site's OMOP database — any
              cohort definable by procedure, diagnosis, demographics or lab value. <code>SELECT</code>-only.
            </li>
            <li>
              <strong>Site evaluation.</strong> Each site's governance lead reviews and can veto,
              independently, through their own Data Access Committee.
            </li>
            <li>
              <strong>Upload the app.</strong> Training code, model config and weights. FLIP scans every
              upload for vulnerabilities before it runs anywhere.
            </li>
            <li>
              <strong>Monitor and collect.</strong> Live per-round metrics per site, then the aggregated
              model with a full provenance trail.
            </li>
          </ol>
          <div class="figure center">
            <img
              :src="deckAsset('ui-cohort-query-to-staging.gif')"
              alt="Animated FLIP walkthrough, recorded end to end: a project with no cohort query yet, an OMOP SQL query written in the cohort query editor, a Run on all trusts action, the per-trust response panel moving each of four trusts from queued to running to a returned count, a live aggregated cohort settling at 5,654 records across four trusts with one count suppressed, then the project page staging each trust and advancing to Awaiting Approval"
              style="max-height: 300px"
            />
            <div class="caption">Steps 1–2: one query, four sites answer, then approval.</div>
          </div>
        </div>
        <p class="small muted" style="margin-top: 0.4em">
          The researcher never accesses site infrastructure and never executes a command on a remote host.
        </p>
      </div>
      <aside class="notes">
        (~45s) This is the map for the demo slides that follow — say it once, then let the screenshots do the
        work. The closing line is the one clinicians and IT directors both care about, from opposite
        directions: no external researcher ever gets a shell on a hospital machine.
      </aside>
    </section>

    <!-- 15 · Demo, step 2 — the approval screen. The one slide that turns "you
         keep control" into something the room can see. -->
    <section>
      <div class="eyebrow">The demo · step 2 · each site decides</div>
      <h2>The approval step, in the product</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1fr 1.1fr">
          <ul class="checklist small">
            <li>The project stays <strong>"Awaiting Approval"</strong> until every named trust decides. Nothing executes locally before that.</li>
            <li>Approval is a <strong>per-trust toggle</strong>, saved explicitly by that site's administrator. One site's switch says nothing about another's.</li>
            <li>Only on approval does FLIP resolve the cohort to DICOM series and cache images from PACS into the site's own XNAT.</li>
            <li>The decision, the decider and the timestamp are written to the audit trail.</li>
          </ul>
          <div class="figure center">
            <img
              :src="deckAsset('ui-approve-project.gif')"
              alt="Animated FLIP project screen: a project awaiting approval, showing the Project Created / Cohort Query / Project Staged / Project Approved stepper, an estimated cohort size across 2 trusts, and a Trust Approval panel where an administrator toggles their own trust from Not Approved and saves"
              style="max-height: 330px"
            />
            <div class="caption">The per-trust approval toggle, in the product.</div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~45s) This slide is where a clinical/governance audience leans in — it's the one screen that
        turns "we promise you keep control" into something they can see. Let the animation loop while
        you talk; the beat that matters is the Trust Approval toggle flipping for one named trust and
        nothing else. If the live demo is running, show the real screen instead.
        Worth saying: approval is also the trigger for image movement — before it, no DICOM has been
        touched. That ordering is the reassurance, not the approval button itself.
      </aside>
    </section>

    <!-- 16 · Demo, steps 3 and 4 — a live federated run with both sites'
         loss curves on one screen. -->
    <section>
      <div class="eyebrow">The demo · steps 3 &amp; 4 · run and watch</div>
      <h2>A federated run, live, per site</h2>
      <div class="slide-body">
        <!-- Rebalanced now the right column carries the QR as well as the figure:
             a wider text column wraps into fewer lines, which is what pays for the
             QR panel's height. -->
        <div class="fig-split" style="--cols: 1.15fr 1.35fr; align-items: center">
          <ul class="checklist small">
            <li>Upload the app bundle, or pick one of the ready-made tutorials.</li>
            <li>Test the setup in a local simulation.</li>
            <li>FLIP schedules the job across every approved site and starts the rounds.</li>
            <li>Live activity feed plus per-site training-loss curves. The researcher watches both ends of the federation from one screen.</li>
            <li>Everything is retained: cohort SQL, site responses, app files, weights, per-round metrics, approvals and logs.</li>
          </ul>
          <!-- Figure and QR share the right-hand grid cell: .fig-split is a
               two-column grid, so a third child would drop the screenshot onto a
               second row. The QR resolves to a read-only snapshot of these same
               Ark+ runs, which is why it sits under the screenshot of them. -->
          <div>
            <div class="figure center">
              <img
                :src="asset('flip-live-training-bdms-kcl.png')"
                alt="FLIP admin UI showing a live federated training run, with the activity feed naming King's College London and Bangkok Dusit Medical Services as the selected trusts, and per-site training-loss curves"
                style="width: 100%; height: auto"
              />
              <div class="caption">Live run: KCL and BDMS training in the same job.</div>
            </div>
            <div class="panel" style="display: flex; align-items: center; gap: 0.8em; margin-top: 0.5em">
              <img
                :src="deckAsset('ark_demo.png')"
                alt="QR code linking to app.flip.aicentre.co.uk/ark_demo, a read-only snapshot of the Ark+ federated experiments in the FLIP interface"
                style="height: 104px; width: 104px; flex: none; border-radius: 6px"
              />
              <div>
                <h3 style="font-size: 0.78em; margin-bottom: 0.1em">Open it yourself</h3>
                <p class="small muted" style="margin: 0; font-size: 0.62em">
                  A read-only snapshot of these Ark+ runs, in the real interface.<br>
                  <code>app.flip.aicentre.co.uk/ark_demo</code>
                </p>
              </div>
            </div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min) Point at the two loss curves and note they belong to machines 9,000 km apart. The
        provenance bullet is the sleeper feature for anyone who has tried to reconstruct what a
        federated experiment actually did six months later — the deployed run keeps the full record,
        which is a different and arguably more useful kind of reproducibility than a released script.
      </aside>
    </section>

    <!-- 17 · Recorded walkthroughs, two workloads side by side. Each clip
         carries its own QR code (both decoded to confirm they point at their own
         video), so the room can take the recordings away even if neither is
         played live. The iframes are sized by their column rather than a fixed
         pixel width, so both stay 16:9 and the pair comes out the same height. -->
    <section>
      <div class="eyebrow">The demo · recorded walkthroughs</div>
      <h2>Complete walkthrough</h2>
      <div class="slide-body">
        <div class="cols" style="--n: 2; gap: 0.8em">
          <div class="figure center" style="width: 100%">
            <iframe
              src="https://www.youtube.com/embed/BH97Yw_lmpw"
              title="FLIP demo — training a chest X-ray classification model"
              style="width: 100%; aspect-ratio: 16 / 9; border: 0; border-radius: 6px"
              allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
              allowfullscreen
            ></iframe>
            <div style="display: flex; align-items: center; justify-content: center; gap: 0.7em; margin-top: 0.45em">
              <img
                :src="deckAsset('demo_video_xray.png')"
                alt="QR code linking to youtube.com/watch?v=BH97Yw_lmpw, the FLIP chest X-ray classification walkthrough"
                style="height: 122px; width: 122px; border-radius: 6px; flex: none"
              />
              <div style="text-align: left">
                <h3 style="font-size: 0.8em; margin-bottom: 0.15em">Chest X-ray classification</h3>
                <div class="caption" style="margin-top: 0">2D, Ark+ fine-tuning<br>youtu.be/BH97Yw_lmpw</div>
              </div>
            </div>
          </div>
          <div class="figure center" style="width: 100%">
            <iframe
              src="https://www.youtube.com/embed/h7fHHkWEuEA"
              title="FLIP demo — training a CT spleen segmentation model"
              style="width: 100%; aspect-ratio: 16 / 9; border: 0; border-radius: 6px"
              allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
              allowfullscreen
            ></iframe>
            <div style="display: flex; align-items: center; justify-content: center; gap: 0.7em; margin-top: 0.45em">
              <img
                :src="deckAsset('demo_video_ct.png')"
                alt="QR code linking to youtube.com/watch?v=h7fHHkWEuEA, the FLIP CT spleen segmentation walkthrough"
                style="height: 122px; width: 122px; border-radius: 6px; flex: none"
              />
              <div style="text-align: left">
                <h3 style="font-size: 0.8em; margin-bottom: 0.15em">CT spleen segmentation</h3>
                <div class="caption" style="margin-top: 0">3D, MONAI bundle<br>youtu.be/h7fHHkWEuEA</div>
              </div>
            </div>
          </div>
        </div>
        <p class="small muted" style="margin-top: 0.4em">
          Two-dimensional classification and three-dimensional organ segmentation share no
          application code — and neither needed a platform change.
        </p>
      </div>
      <aside class="notes">
        (~1 min) The reuse claim made concrete: 3D organ segmentation shares no application code with
        2D classification, and neither required a platform change. The spleen clip is filmed from the
        site administrator's side rather than the researcher's, which shows the governance half of the
        workflow. If time is tight, play neither and point at the QR codes — the room can watch both
        later, and the job-types slide has already made the point.
      </aside>
    </section>

    <!-- 18 · The team — engineering, clinical and governance people in one
         group, which is itself part of why the deployment happened. -->
    <section>
      <div class="eyebrow">Building FLIP together</div>
      <h2>Our team</h2>
      <div class="slide-body">
        <div class="cols" style="display: block; width: fit-content; margin: 0 auto">
          <img
            :src="sharedAsset('shared/seb_kcl.png')"
            alt="Professor Sébastien Ourselin, Head of School, School of Biomedical Engineering & Imaging Sciences"
            style="max-height: 165px"
          />
        </div>
        <div class="fig-split" style="--cols: 4fr 3fr; margin-top: 0.3em; gap: 1.2em; align-items: start">
          <div>
            <div class="eyebrow" style="text-align: center; color: #eb2f2d; font-weight: 600; font-size: large;">King's College London</div>
            <div style="display: flex; justify-content: center; flex-wrap: wrap; gap: 0.4em; margin-top: 0.3em">
              <img :src="sharedAsset('shared/jorge.png')" alt="Dr M. Jorge Cardoso, Group Lead & Reader" style="height: 170px" />
              <img :src="sharedAsset('shared/rafa_kcl.png')" alt="Rafael Dias, Senior AI Engineer on Foundational Models for Healthcare" style="height: 170px" />
              <img :src="sharedAsset('shared/alex_kcl.png')" alt="Alexandre Triay Bagur, Senior AI Engineer" style="height: 170px" />
              <img :src="sharedAsset('shared/virginia_kcl.png')" alt="Virginia Fernandez, Research Associate" style="height: 170px" />
            </div>
          </div>
          <div>
            <div class="eyebrow" style="text-align: center; color: #005EB8; font-weight: 600; font-size: large;">Guy's and St Thomas' Trust</div>
            <div style="display: flex; justify-content: center; flex-wrap: wrap; gap: 0.4em; margin-top: 0.3em">
              <img :src="sharedAsset('shared/joe_gstt.png')" alt="Joe Zhang, Head of Data Science" style="height: 170px; border-radius: 8px" />
              <img :src="sharedAsset('shared/lawrence_gstt.png')" alt="Lawrence Adams, Lead Analytics Engineer" style="height: 170px; border-radius: 8px" />
              <img :src="sharedAsset('shared/martin_gstt.png')" alt="Martin Chapman, Lead NLP Engineer" style="height: 170px; border-radius: 8px" />
            </div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~20s) Fast slide — put faces to the platform and name the collaboration. Worth one sentence
        that this is engineering, clinical and governance people in the same team, which is itself part
        of why the deployment happened.
      </aside>
    </section>

    <!-- 19 · Close. Stays up through Q&A so the three QR codes remain on
         screen. -->
    <section class="title-slide center">
      <div class="eyebrow">Federated Learning Interoperability Platform</div>
      <h1>Send the model to the data</h1>
      <div class="slide-body">
        <p class="subtitle">
          Open source, Apache 2.0, running today across two continents. Come and use it.
        </p>
        <div class="cols" style="--n: 3; gap: 0.7em; max-width: 24em; margin: 0.8em auto 0; align-items: stretch">
          <div>
            <p class="figure-placeholder__desc" style="margin-bottom: 0.35em">Contact me</p>
            <div class="figure-placeholder compact" style="aspect-ratio: 1 / 1; min-height: 0; padding: 0.55em">
              <img
                :src="pydataAsset('me.png')"
                alt="https://garciadias.github.io/#/"
                title="Contact me"
                style="display: block; width: 100%; height: 100%; object-fit: cover; border-radius: 6px"
              />
            </div>
          </div>
          <div>
            <p class="figure-placeholder__desc" style="margin-bottom: 0.35em">FLIP repository</p>
            <div class="figure-placeholder compact" style="aspect-ratio: 1 / 1; min-height: 0; padding: 0.55em">
              <img
                :src="pydataAsset('FLIP_REPO.png')"
                alt="https://github.com/londonaicentre/FLIP"
                title="FLIP repository"
                style="display: block; width: 100%; height: 100%; object-fit: cover; border-radius: 6px"
              />
            </div>
          </div>
          <div>
            <p class="figure-placeholder__desc" style="margin-bottom: 0.35em">AI Centre team</p>
            <div class="figure-placeholder compact" style="aspect-ratio: 1 / 1; min-height: 0; padding: 0.55em">
              <img
                :src="pydataAsset('ai_centre.png')"
                alt="https://www.aicentre.co.uk/about/meet-the-team"
                title="AI Centre team"
                style="display: block; width: 100%; height: 100%; object-fit: cover; border-radius: 6px"
              />
            </div>
          </div>
        </div>
        <p class="byline">Questions — especially the awkward ones.</p>
      </div>
      <aside class="notes">
        (~20s, then 5 min Q&amp;A) Land on the one-line identity: open source, Apache 2.0, no lock-in,
        already running. Leave this slide up through Q&amp;A so the QR codes stay on screen.
        Likely questions to have ready: differential privacy and formal guarantees; what happens when
        the funding ends; how FLIP compares to Flower's own deployment tooling; whether a site can
        leave; how long onboarding really takes; and why synthetic data.
      </aside>
    </section>
  </RevealDeck>
</template>
