<script setup>
import DeckBrand from '@/components/DeckBrand.vue';
import RevealDeck from '@/components/RevealDeck.vue';

// 25-minute talk for "From Zero to Hero: Federated AI in Healthcare Systems"
// (St John's College, Cambridge). Audience: clinicians, researchers and people
// new to federated learning — so the argument is deliberately platform-first
// rather than algorithm-first.
//
// Narrative, in four movements (slide numbers in the comments below match the
// order of the sections in this file, 0-20):
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
//   18-20 what's next: the platforms other groups are building and the offer to
//         converge with them, the team, and the close.
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
        <p class="subtitle">
          Federated Learning & Interoperability Platform
        </p>
        <p class="byline">
          Rafael Garcia-Dias
        </p>
        <p class="venue-note">github.com/londonaicentre/FLIP · Apache 2.0</p>
      </div>
      <aside class="notes">
        (~15s) (at 0:00) 
      </aside>
    </section>

    <!-- 1 · Roadmap. Four beats, so the room knows the demo is coming and can
         hold its questions until the end. -->
    <section>
      <div class="eyebrow">FLIP: Federated Learning Interoperability Platform</div>
      <h2>What you'll get in the next ~20 minutes</h2>
      <div class="slide-body medium">
        <ol class="contribs">
          <li><strong>The gap</strong> Why are the FL frameworks ready, but not yet mainstream in hospitals?</li>
          <li><strong>The platform</strong> What FLIP is, what it's built on, and what it deliberately doesn't do?</li>
          <li><strong>The demo</strong> A cross-continental deployment, cohort query, governance, training run, results.</li>
          <li><strong>What's next</strong> How to join the community and contribute?</li>
        </ol>
      </div>
      <aside class="notes">
        (~45s) (at 0:15)
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
        (~45s) (at 1:00) Why FLIP focus on imaging?
      </aside>
    </section>

    <!-- 3 · The problem, 2 of 2  -->
    <section>
      <div class="eyebrow">The problem · 2 of 2 · the thesis of this talk</div>
      <h2>On the shoulders of giants; Connecting institutions.</h2>
      <div class="slide-body">
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
        (~1 min) (at 1:45) Why FLIP exists?
      </aside>
    </section>

    <!-- 4 · Architecture, 1 of 2  -->
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
        (~1 min) (at 2:45) Central hub and Nodes
      </aside>
    </section>

    <!-- 5 · Architecture, 2 of 2  -->
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
        (~3 min) (at 3:45) UI, CH, Trust, APIs, FL clients and servers
      </aside>
    </section>

    <!-- 6 · Microservice inventory. -->
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
        (~15s) (at 6:45) For people to see offline, don't explain much
      </aside>
    </section>

    <!-- 9 · Deployment modes — three objections (how, what hardware, what network) answered before anyone raises them. -->
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
        (~1min) (at 7:00) How easy is to deploy a new node?
      </aside>
    </section>

    <!-- 7 · The outbound-only mechanism -->
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
            <div class="stat-label">API verbs; the hub's entire request vocabulary</div>
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
        (~1 min) (at 8:00) Simplify it for IT and Trust security teams
      </aside>
    </section>

    <!-- 8 · Where the guardrails physically execute. -->
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
        (~1 min) (at 9:00) More security guarantees
      </aside>
    </section>

    <!-- 10 · Governance as architecture -->
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
        (~1 min) (at 10:00) Governance is architecture.
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
        (~45s) (at 11:00) Focus on the researcher's view. The platform needs to be easy to use.
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
        (~5s) (at 11:45)
      </aside>
    </section>

    <!-- 13 · The deployment being demoed: the estate at both ends, and the task
         it was pointed at, stacked in one figure card. -->
    <section>
      <div class="eyebrow">The demo · the deployment behind it</div>
      <h2>Cross-continental deployment</h2>
      <div class="center">
        <p class="small">
          KCL and BDMS (Bangkok Dusit Medical Services)
        </p>
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
        (~45s) (at 11:50) BDMS example. 

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
        (~45s) (at 12:35) This is the map for the demo slides that follow — say it once, then let the screenshots do the
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
          <ul class="dotlist small">
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
        (~45s) (at 13:20) 
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
        <div class="fig-split" style="--cols: 1.05fr 1.3fr; align-items: center">
          <div class="contribs small">
            <ul class="dotlist">
              <li>Upload the app bundle, or pick one of the ready-made tutorials.</li>
              <li>Test the setup in a local simulation.</li>
              <li>FLIP schedules the job across every approved site and starts the rounds.</li>
              <li>Live activity feed plus per-site training-loss curves. The researcher watches both ends of the federation from one screen.</li>
              <li>Everything is retained: cohort SQL, site responses, app files, weights, per-round metrics, approvals and logs.</li>
            </ul>
          </div>
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
        (~1 min) (at 14:05) 
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
        (~1 min) (at 15:05) 
      </aside>
    </section>

    <!-- 18 · The wider ecosystem, opening the "what's next" movement promised on
         slide 1. Deliberately a roll-call and not the capability matrix it was
         drawn from: the point is that these exist and were built for different
         jobs, so a per-platform scorecard would invite exactly the line-by-line
         argument this slide is trying to avoid. Grouped four ways, which is the
         only claim being made about them. Three columns rather than four —
         peer-to-peer (2 names) and analytics (3) stack into the middle one, so
         the widest group (the engines, 8) gets a column to wrap in and the
         panels come out roughly level. One short line per group and nothing more:
         at a third of the slide's width, a sentence of any length wraps to four
         or five lines and the panels grow to pay for it. The FLIP pill is marked
         .ours so the room can find us in the list without the slide having to
         rank anything.

         The .here dots are the cooperation argument, and every one was checked
         against this workshop's own programme (FromZeroToHero_programme_v4.pdf)
         and the platforms' papers:
           Flower    — Nicholas Lane (co-founder/CSO, chairs panel 2) and Javier
                       Fernandez-Marques (Flower SuperGrid demo, tutorial 2);
                       Flower Labs is also the exclusive sponsor.
           PharosAI  — Gregory Verghese (co-founder/CTO, chairs panel 1), Concetta
                       Piazzese and Bangxiang Guan (talks), Olaoluwa Ademolu (on
                       panel 2 with us). KCL + QMUL + GSTT + Barts, so partly the
                       same institutions as the team slide. Not in the comparison
                       table this slide was drawn from — it belongs there.
           FLA³      — Fan Zhang (lead organiser, first author) and Michael Roberts
                       (BloodCounts! lead); arXiv 2603.10063 also lists Kreuter,
                       Fernandez-Marques, Verghese, Lane and Schönlieb, which is
                       seven of the nine organisers.
           vantage6  — André Dekker (keynote 1) leads Clinical Data Science
                       Maastricht, which co-developed vantage6 and originated the
                       Personal Health Train concept behind it.
         Carsten Maple, Mark Durkee, John Smith, Shiqiang Wang, Bing Luo and Huw
         Day work on federated security, government pilots and FL research rather
         than a named platform, so they are deliberately not dotted. -->
    <section>
      <div class="eyebrow">What's next · the wider ecosystem</div>
      <h2>An ecosystem, not a race</h2>
      <div class="slide-body">
        <p class="small">
          FLIP is one of many federated learning platforms. Each was built for a different purpose, and each has its own
          strengths. The point is not to compare them, but to show that the ecosystem is rich and that agreeing on
          interfaces is cheaper than rebuilding platforms.
        </p>
        <div class="cols compact" style="--n: 3; gap: 0.8em; align-items: start">
          <div class="panel flip">
            <h3>Deployed platforms</h3>
            <p class="small muted">End to end, inside hospital estates. Our closest neighbours.</p>
            <div>
              <span class="pill ours here">FLIP</span><span class="pill here">PharosAI</span
              ><span class="pill here">FLA<sup>3</sup></span><span class="pill">CODA</span
              ><span class="pill">FeatureCloud</span><span class="pill">JIP / Kaapana</span
              ><span class="pill">MedPerf</span>
            </div>
          </div>
          <div>
            <div class="panel">
              <h3>Federated analytics</h3>
              <p class="small muted">Statistics across sites, no imaging path.</p>
              <div>
                <span class="pill here">vantage6</span><span class="pill">TrainTracks</span
                ><span class="pill">DataSHIELD</span>
              </div>
            </div>
            <div class="panel" style="margin-top: 0.55em">
              <h3>Peer-to-peer</h3>
              <p class="small muted">No central coordinator. The one thing FLIP does not do.</p>
              <div>
                <span class="pill">Swarm Learning</span><span class="pill">Fedstellar</span>
              </div>
            </div>
          </div>
          <div class="panel">
            <h3>FL engines</h3>
            <p class="small muted">The training core. FLIP runs on two of these and competes with none.</p>
            <div>
              <span class="pill here">Flower</span><span class="pill">NVIDIA FLARE</span
              ><span class="pill">OpenFL</span><span class="pill">Fed-BioMed</span
              ><span class="pill">Substra</span><span class="pill">FEDn</span
              ><span class="pill">FATE</span><span class="pill">FedML</span>
            </div>
          </div>
        </div>
        <p class="small muted" style="margin-top: 0.5em">
          <span class="here-dot"></span>someone building it is in this room this week
        </p>
      </div>
      <aside class="notes">
        (~1 min) (at 16:05) Generous tone, and no scorecard: read the room's own work back to it — four
        different jobs, not four attempts at the same one. Concede the peer-to-peer column
        outright: FLIP has a central hub and always will, so that group is solving something we
        are not. Concede the analytics column's governance maturity too.

        Then work the dots, because they are the whole point of the slide and they are all in the
        room: Flower — Nic Lane, chairing our own panel, and Javier Fernandez-Marques, who demos
        SuperGrid and runs tutorial 2. PharosAI — Gregory Verghese chairing panel 1, plus Concetta
        Piazzese, Bangxiang Guan, and Olaoluwa Ademolu sitting on panel 2 next to me; and it is
        KCL, Guy's and Barts, so partly the same institutions as our own team slide. FLA³ — Fan
        Zhang and Michael Roberts, who organised this workshop. vantage6 — André Dekker's
        department co-developed it, and this morning's keynote is his.

        The strongest single fact, worth saying out loud: seven of the nine people who organised
        this workshop are co-authors on FLA³, one of the platforms on this slide — Flower Labs and
        PharosAI among them. Nobody in this field is actually competing; we just have not agreed
        the seams yet. That is the invitation.

        If asked for the detailed comparison, say it is written up as a capability matrix per
        layer and offer it afterwards rather than putting it on screen.
      </aside>
    </section>

    <!-- 19 · The team — engineering, clinical and governance people in one
         group, which is itself part of why the deployment happened. -->
    <section>
      <div class="eyebrow">Building FLIP together</div>
      <h2>FLIP team</h2>
      <div class="slide-body">
        <!-- Every headshot sits in a .photo-chip, which is what rounds its corners:
             these PNGs are flat cards with the name, role and a background colour
             baked in, so the softening has to be clipped on by the wrapper. -->
        <div class="people" style="margin-top: 0">
          <div class="photo-chip">
            <img
              :src="sharedAsset('shared/seb_kcl.png')"
              alt="Professor Sébastien Ourselin, Head of School, School of Biomedical Engineering & Imaging Sciences"
              style="max-height: 165px"
            />
          </div>
        </div>
        <div class="fig-split" style="--cols: 4fr 3fr; margin-top: 0.3em; gap: 1.2em; align-items: start">
          <div>
            <div class="eyebrow" style="text-align: center; color: #eb2f2d; font-weight: 600; font-size: large;">King's College London</div>
            <div class="people">
              <div class="photo-chip">
                <img :src="sharedAsset('shared/jorge.png')" alt="Dr M. Jorge Cardoso, Group Lead & Reader" style="height: 170px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/rafa_kcl.png')" alt="Rafael Dias, Senior AI Engineer on Foundational Models for Healthcare" style="height: 170px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/alex_kcl.png')" alt="Alexandre Triay Bagur, Senior AI Engineer" style="height: 170px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/virginia_kcl.png')" alt="Virginia Fernandez, Research Associate" style="height: 170px" />
              </div>
            </div>
          </div>
          <div>
            <div class="eyebrow" style="text-align: center; color: #005EB8; font-weight: 600; font-size: large;">Guy's and St Thomas' Trust</div>
            <div class="people">
              <div class="photo-chip">
                <img :src="sharedAsset('shared/joe_gstt.png')" alt="Joe Zhang, Head of Data Science" style="height: 170px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/lawrence_gstt.png')" alt="Lawrence Adams, Lead Analytics Engineer" style="height: 170px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/martin_gstt.png')" alt="Martin Chapman, Lead NLP Engineer" style="height: 170px" />
              </div>
            </div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~20s) (at 17:05) 
      </aside>
    </section>

    <!-- 20 · Close. Stays up through Q&A so the three QR codes remain on
         screen. -->
    <section class="title-slide center">
      <div class="eyebrow">Federated Learning Interoperability Platform</div>
      <h1>Send the model to the data</h1>
      <div class="slide-body">
        <p class="subtitle">
          Open source, Apache 2.0. Use it, change it, contribute.
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
      </div>
      <aside class="notes">
        (~20s, then 5 min Q&amp;A) (at 17:25) 
      </aside>
    </section>
  </RevealDeck>
</template>
