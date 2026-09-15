<script setup>
import DeckBrand from '@/components/DeckBrand.vue';
import RevealDeck from '@/components/RevealDeck.vue';

// 20-minute talk for NVIDIA FLARE Day 2026 (16 September, US + EMEA main
// event, 10:10-10:30 PST). Derived from FlipReusableFlPlatformDeck.vue (the
// Cambridge workshop talk) but retargeted in four ways:
//
//   1. Audience. Cambridge was clinicians and researchers new to FL, so that
//      deck had to argue that federated learning matters at all. This room is
//      FLARE practitioners — they are already convinced. The Gboard privacy
//      slide and the "why imaging" framing are cut; the deck opens on the gap
//      between a framework and a hospital, which is the thing this audience
//      has personally hit.
//   2. FLARE is the host. NVIDIA organises this event, so the FLARE
//      relationship gets its own slide (§4) rather than a passing mention:
//      what FLIP delegates to FLARE, what it adds, and the MONAI bundle path.
//   3. The demo stops being a separate act. Cambridge had a "Part three ·
//      Demo" divider slide and a self-contained ~7 minute demo section,
//      because the organisers asked for it. Here the product screens are woven
//      into the platform argument as evidence for the claim being made at that
//      moment — the cohort-query GIF sits on the cohort slide, the approval
//      GIF on the governance slide. No gear change, no divider.
//   4. The ecosystem slide is rebuilt around this event. The Cambridge version
//      dotted people who were in the room at St John's. That roll-call is
//      meaningless here, so it is fused with a short cut of the capability
//      comparison (PlatformComparisonDeck.vue) and re-dotted against the
//      FLARE Day programme — every name checked against the published agenda
//      on events.nvidia.com/flare-day-2026 (see §11 comment for the mapping).
//
// Four slides were later folded in from ".data/2026-09-10 NECTEC.odp" — an
// earlier FLIP deck whose security (§13), MONAI Application Package (§15),
// FLIP-SIDE (§14) and AWS Landing Zone (§6) slides carried material this talk
// was missing. The compliance-framework names, the six sovereignty guarantees,
// the Ark+ run specifics (50 FedAvg rounds, ~27 KB/client/round, 6,885 trainable
// params) and the AWS architecture diagram all come from there; everything else
// was rebuilt as native markup so it matches this deck's styling.
//
// Slide budget for a 20-minute slot: 16 slides, ~17 min of speaking, ~3 min Q&A.
//
// deckAsset()   → this deck's folder: the intro globe (copied so the deck is
//                 self-contained if the Cambridge folder is ever pruned).
// reusableAsset() → the Cambridge deck's folder: bridge artwork, UI
//                 recordings, QR codes. Reused rather than duplicated.
// asset()       → flip-maturity-pitch-2026 image folder.
// sharedAsset() → public/presentations/ (team headshots).
// pydataAsset() → QR codes built for the PyData lightning talk.
const deckAsset = (name) => `${import.meta.env.BASE_URL}presentations/flare-day-2026/${name}`
const reusableAsset = (name) => `${import.meta.env.BASE_URL}presentations/flip-reusable-fl-platform/${name}`
const asset = (name) => `${import.meta.env.BASE_URL}presentations/flip-maturity-pitch-2026/${name}`
const sharedAsset = (name) => `${import.meta.env.BASE_URL}presentations/${name}`
const pydataAsset = (name) => `${import.meta.env.BASE_URL}presentations/pydata-london-flip-lightning-2026/${name}`
</script>

<template>
  <RevealDeck
    :options="{ center: true }"
    :style="{ '--hero-poster': `url(${deckAsset('intro-globe-poster.jpg')})` }"
  >
    <template #chrome>
      <DeckBrand :qr="reusableAsset('presentation.png')" />
    </template>

    <!-- 1 · Title. Same globe hero as Cambridge; the eyebrow and subtitle carry
         the registered session title and the three partner institutions, which
         is how the programme lists this talk. -->
    <section
      class="title-slide center video-hero"
      :data-background-video="`${deckAsset('intro-globe.mp4')},${deckAsset('intro-globe.webm')}`"
      data-background-video-loop
      data-background-video-muted
      data-background-size="cover"
      data-background-color="#ffffff"
    >
      <div class="eyebrow">NVIDIA FLARE Day 2026 · 16 September · US + EMEA main event</div>
      <h1>FLIP: an open-source federated learning platform for healthcare</h1>
      <div class="slide-body">
        <p class="subtitle">
          From multi-institutional research to real NHS deployment
        </p>
        <p class="byline">
          Rafael Garcia-Dias · Senior AI Engineer, King's College London
        </p>
        <p class="venue-note">
          London AI Centre · King's College London · Guy's and St Thomas' NHS Foundation Trust<br>
          github.com/londonaicentre/FLIP · Apache 2.0
        </p>
      </div>
      <aside class="notes">
        (~20s) (at 0:00) Name the three institutions — this audience does not know
        the London AI Centre, and the NHS trust is what makes the deployment claim
        credible. Then go straight to the gap: you have 20 minutes and this room
        already believes in FL.
      </aside>
    </section>

    <!-- 2 · The gap. Rewritten for this audience: Cambridge opened by arguing
         that medical imaging data is valuable and stuck. A FLARE Day room does
         not need that argument — they need the specific claim that the
         framework is not the missing piece. The NHS imaging chart is kept
         because it sizes the prize in one glance, but the text around it is
         about the platform gap, not about why FL exists. -->
    <section>
      <div class="eyebrow">The gap · why a framework is not yet a deployment</div>
      <h2>FLARE solved the hard part. Then the hospital said no.</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1.3fr 0.95fr; margin-top: 0.2em; align-items: start">
          <div>
            <p class="small">
              Every major NHS trust has a PACS archive going back decades. The federated core to train
              on it has existed for years. The blocker was never the aggregation strategy.
            </p>
            <ul class="crosslist small">
              <li><strong>The network team, not the data scientist</strong> An inbound port is a governance
                conversation that outlasts the project.</li>
              <li><strong>The cohort does not exist yet</strong> "Patients with a chest X-ray and a
                positive culture" is a query against each site's own schema, before any training starts.</li>
              <li><strong>Approval is per site, per project</strong> And it has to gate execution, not
                just the UI.</li>
            </ul>
            <p class="small" style="margin-top: 0.5em">
              Everything in that list sits <em>between</em> the hospital and FLARE — and every
              collaboration rebuilds it by hand.
            </p>
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
        (~1 min) (at 0:20) The one-liner that lands in this room: we did not have a
        federated learning problem, we had an everything-around-federated-learning
        problem. Do not relitigate why FL matters — they run FLARE already.
      </aside>
    </section>

    <!-- 3 · The thesis, with the bridge artwork. Kept from Cambridge because it
         is the single clearest statement of what FLIP is, but the panel now
         names FLARE first since FLARE is the host framework here. -->
    <section>
      <div class="eyebrow">The thesis of this talk</div>
      <h2>On the shoulders of giants: connecting institutions</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1.65fr 1fr; margin-top: 0; align-items: start">
          <div>
            <div class="panel flip">
              <h3>The federated core is ready</h3>
              <p class="small">
                <strong>NVIDIA FLARE</strong> and <strong>Flower</strong> made it something you can
                build and depend on. Aggregation strategies, client/server orchestration, secure
                transport, real communities behind both.
              </p>
            </div>
            <p class="small" style="margin-top: 0.5em">
              FLIP is the bridge <em>between a hospital and the FL network</em>
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
              :src="reusableAsset('flip_bridge.png')"
              alt="Illustration: two hospitals on opposite sides of a river valley, each with its own PACS, EHR and governance office. A stone bridge labelled FLIP — Federated Learning Interoperability Platform — spans the gap between them, carried on two piers marked Flower and NVIDIA FLARE"
              style="max-width: 100%; max-height: 470px"
            />
          </div>
        </div>
      </div>
      <aside class="notes">
        (~45s) (at 1:20) Point at the piers in the drawing — FLARE and Flower hold
        the bridge up. We are not competing with anything in this event's stack;
        we are the span that sits on top of it.
      </aside>
    </section>

    <!-- 4 · NEW · The FLARE + MONAI slide. Requested explicitly, and the one
         slide most specific to this audience. Laid out as three columns of equal
         height rather than a two-column split with a stacked right side, which
         ran 150px past the 720px box. The MONAI panel covers both directions:
         bundles in (the 3D segmentation job type runs a stock MONAI bundle) and
         packages out (trained models wrapped as MONAI Application Packages for
         the reading environment) — the second half comes from the NECTEC deck's
         §15 and closes the loop the deck would otherwise leave open. -->
    <section>
      <div class="eyebrow">Where FLIP sits on the NVIDIA stack</div>
      <h2>FLARE for orchestration, MONAI for the medical imaging</h2>
      <div class="slide-body">
        <div class="cols" style="--n: 3; gap: 0.6em; align-items: start">
          <div class="panel flip">
            <h3>Delegated to FLARE</h3>
            <ul class="dotlist small" style="margin-top: 0.15em; margin-bottom: 0">
              <li>Client/server orchestration, round control</li>
              <li>Aggregation strategies (FedAvg and beyond)</li>
              <li>mTLS transport, startup-kit provisioning</li>
              <li>Event hooks our audit records hang off</li>
            </ul>
          </div>
          <div class="panel">
            <h3>Added by FLIP</h3>
            <ul class="dotlist small" style="margin-top: 0.15em; margin-bottom: 0">
              <li>Cohort query on each site's OMOP database</li>
              <li>PACS retrieval into site-local XNAT</li>
              <li>Per-site approval that gates execution</li>
              <li>Outbound-only transport and cert provisioning</li>
            </ul>
          </div>
          <div class="panel">
            <h3>MONAI, both ways</h3>
            <ul class="dotlist small" style="margin-top: 0.15em; margin-bottom: 0">
              <li><strong>In:</strong> the segmentation job type runs an
                <strong>unmodified MONAI bundle</strong>, own transforms and config</li>
              <li><strong>Out:</strong> trained models leave as <strong>MONAI Application
                Packages</strong> for the reading environment</li>
            </ul>
          </div>
        </div>
        <p class="small muted" style="margin: 0.4em 0 0; text-align: center">
          FLIP provisions <strong>one isolated FLARE network per run</strong>, one client per Trust —
          and resolves the cohort to on-site DICOM before the first round.
        </p>
      </div>
      <aside class="notes">
        (~1 min 30) (at 2:05) The slide this room came for. Three beats. One — we
        did not fork FLARE and we did not wrap it in something that hides it; the
        provisioning model and the filter hooks are used as designed, and our audit
        records hang off FLARE events. Two — everything in the middle column is
        the stuff FLARE quite reasonably does not do, because it is hospital
        integration, not federated learning. Three — MONAI both ways: the
        segmentation job trains a stock bundle with no FLIP-specific code in it,
        and a trained model leaves as a MONAI Application Package so it lands back
        in the reading environment where a radiologist actually works. If someone
        asks about MONAI FL / the FLARE-MONAI integration, the honest answer is
        that we use the bundle format rather than MONAI's own FL client, because
        the client side has to go through our trust-node data path.
      </aside>
    </section>

    <!-- 5 · Architecture. The detailed diagram from Cambridge. This is the
         slide that answers "what did you actually build", and for this audience
         the interesting part is that only two links cross the node boundary. -->
    <section>
      <div class="eyebrow">The platform · how it physically works</div>
      <h2>One Central Hub, one Trust Node per site</h2>
      <div class="slide-body">
        <div class="figure center" style="width: 100%">
          <img
            :src="reusableAsset('flip_architecture-flip_architecture.png')"
            alt="Detailed FLIP architecture: a cloud-hosted FLIP Hub holding the FLIP UI, Hub API, FL API and FL Server; on the right, a hospital's EHR and PACS feeding a structured OMOP database and an XNAT imaging archive inside a secure enclave, each read by its own Data API and Imaging API, with the Trust API orchestrating both. Only two links cross the FLIP node boundary: Hub API to Trust API, and FL Server to FL Client"
            style="max-height: 500px; max-width: 100%"
          />
        </div>
        <p class="small muted" style="margin-top: 0.25em; text-align: center; margin-bottom: 0">
          Only two links cross the trust boundary: <code>Hub API → Trust API</code>, and
          <code>FL Server → FL Client</code>. Both opened from inside.
        </p>
      </div>
      <aside class="notes">
        (~1 min 15) (at 3:35) Walk it once, left to right, then land on the two
        crossings. Everything else on this diagram is inside the hospital's own
        enclave and stays there. Note the FL Server / FL Client pair is FLARE's —
        we provision it per run.
      </aside>
    </section>

    <!-- 6 · The outbound-only mechanism. The single most transferable piece of
         engineering in the talk for a FLARE audience: how you get a FLARE
         client running inside an NHS trust without an inbound firewall rule.
         The session abstract promises "certificate provisioning" as a lesson
         learned, so the key-handling bullet is expanded relative to Cambridge. -->
    <section>
      <div class="eyebrow">The platform · what we asked the network team for</div>
      <h2>Outbound-only, or the answer is no</h2>
      <div class="slide-body">
        <div class="kpis" style="--n: 4">
          <div class="panel center">
            <div class="stat" style="font-size: 1.5em">0</div>
            <div class="stat-label">inbound ports on the trust host, and no port 22 open either</div>
          </div>
          <div class="panel center">
            <div class="stat" style="font-size: 1.5em">5s</div>
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
            <li>The FL clients hold <strong>no Central Hub credentials at all</strong>. Compromising a client does not give an attacker access to the hub.</li>
            <li><strong>Certificate provisioning</strong> Registration mints the site's startup kit and keys; the hub keeps only a <strong>SHA&#8209;256 hash</strong> of the API key, and the plaintext never leaves the site's own kit file.</li>
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
        (~1 min 15) (at 4:50) This is the slide people will ask about afterwards.
        The framing that works with NHS network teams: we are not asking you to
        trust us, we are asking you to let a machine you own make an outbound HTTPS
        call to a published endpoint, exactly like every other thing on your
        network. The 6-verb vocabulary is the security argument — the whole
        attack surface fits on one line of a slide.
      </aside>
    </section>

    <!-- 7 · Deployment modes. Answers "will it fit our stack" before it is
         asked. This audience includes people who deploy FLARE across very
         heterogeneous estates (Mayo's talk earlier in the morning makes exactly
         this point), so the hardware floor is the number worth saying. -->
    <section>
      <div class="eyebrow">The platform · will it fit our stack?</div>
      <h2>Low-barrier deployment on Trust Nodes</h2>
      <div class="slide-body">
        <div class="cols" style="--n: 3">
          <div class="panel">
            <h3>Deploy how you already operate</h3>
            <p class="small">
              Microservices deployable via <strong>Docker Compose</strong>, a <strong>Kubernetes</strong>
              Helm chart, or <strong>cloud</strong> infrastructure-as-code. A site picks the one its IT
              team already runs.
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
        <div class="flip-note">
          <span class="flip-tag">Lesson learned</span> The binding constraint on joining a federation
          was never GPU capacity. It was whether the site's IT team recognised the deployment shape.
        </div>
      </div>
      <aside class="notes">
        (~45s) (at 6:05) Charlie Qin from Mayo Clinic Platform made this exact point
        earlier this morning — deployment across highly variable site environments,
        from mature cloud teams to resource-constrained local IT. Worth naming him
        if the timing works; it is a genuine convergence, not a courtesy.
      </aside>
    </section>

    <!-- 8 · Cohort definition, WITH the product. First of the three woven demo
         slides. Cambridge had a "the researcher's four steps" slide and then a
         separate demo section; here the recording of steps 1-2 sits directly on
         the slide making the claim, so the screen capture is evidence rather
         than a change of act. The four-step list is compressed into the left
         column because the GIF carries the first two steps visually. -->
    <section>
      <div class="eyebrow">The platform · step one, in the product</div>
      <h2>The cohort is a query, not a data request</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1fr 1.15fr; align-items: center">
          <div>
            <p class="small">
              One SQL query runs against <strong>each site's own OMOP database</strong>. Any cohort
              definable by procedure, diagnosis, demographics or lab value. The researcher gets
              <strong>counts back, never rows</strong>.
            </p>
            <ul class="dotlist small">
              <li>The authoritative validator is the <strong>trust's</strong>: single statement,
                <code>SELECT</code> only, pinned to the <code>omop</code> schema, re-emitted from the
                checked syntax tree and run as a <strong>read-only role</strong>.</li>
              <li>Each site sets a <strong>disclosure floor</strong> (default 10). Below it, the
                refusal is identical to an empty result — so narrowing a query until one patient
                matches tells you nothing.</li>
              <li>No keyword denylist anywhere. We parse.</li>
            </ul>
          </div>
          <div class="figure center">
            <img
              :src="reusableAsset('ui-cohort-query-to-staging.gif')"
              alt="Animated FLIP walkthrough, recorded end to end: a project with no cohort query yet, an OMOP SQL query written in the cohort query editor, a Run on all trusts action, the per-trust response panel moving each of four trusts from queued to running to a returned count, a live aggregated cohort settling at 5,654 records across four trusts with one count suppressed, then the project page staging each trust and advancing to Awaiting Approval"
              style="max-height: 340px"
            />
            <div class="caption">One query, four sites answer, one count suppressed by the floor.</div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min 15) (at 6:50) Let the recording run while you talk — it is a real
        capture, not a mock. Call out the suppressed count when it appears: that is
        the disclosure floor doing its job on screen, and it is the most concrete
        privacy guarantee in the talk. The parse-don't-denylist line usually gets a
        nod from the security people in the room.
      </aside>
    </section>

    <!-- 9 · Security and sovereignty, WITH the approval screen. Second woven
         demo slide. Expanded from the NECTEC deck's §13 ("Security and
         sovereignty of NHS data"), which carried six numbered guarantees and
         the named compliance frameworks — strictly richer than the Cambridge
         version this replaces, so the two are merged here: the NECTEC
         guarantees carry the argument and the recorded per-trust approval
         toggle is the evidence for guarantee 01.
         Compliance names are kept verbatim from that deck. They are specific
         and checkable, and this is an NHS-regulatory audience. -->
    <section>
      <div class="eyebrow">The platform · step two, in the product</div>
      <h2>Security and sovereignty of NHS data</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1.25fr 1fr; align-items: center">
          <div>
            <p class="small" style="margin-bottom: 0.25em">
              No project can use a Trust's data until that Trust approves it — <strong>the data and
              the approval both live on the Trust's own infrastructure</strong>.
            </p>
            <div class="cols compact" style="--n: 2; gap: 0.5em">
              <ul class="dotlist small" style="margin: 0">
                <li><strong>Per-project veto.</strong> Approved separately per Trust. Declining one
                  affects nothing else; the hub cannot bypass the gate.</li>
                <li><strong>Data residency.</strong> Only aggregate statistics, filtered model updates
                  and task status leave the Trust.</li>
                <li><strong>Structural guarantees.</strong> Outbound polling only — the hub has no
                  route in, so there is no path to trust.</li>
              </ul>
              <ul class="dotlist small" style="margin: 0">
                <li><strong>National Data Opt-Out.</strong> Applied where data is prepared, so it
                  propagates to every query and model.</li>
                <li><strong>Access governance.</strong> Default-deny roles, MFA on every request,
                  per-project membership.</li>
                <li><strong>Ethical research.</strong> HRA REC-approved Research Database, overseen
                  by a Data Allocation Committee.</li>
              </ul>
            </div>
          </div>
          <div>
            <div class="figure center">
              <img
                :src="reusableAsset('ui-approve-project.gif')"
                alt="Animated FLIP project screen: a project awaiting approval, showing the Project Created / Cohort Query / Project Staged / Project Approved stepper, an estimated cohort size across 2 trusts, and a Trust Approval panel where an administrator toggles their own trust from Not Approved and saves"
                style="max-height: 205px"
              />
              <div class="caption">Guarantee 01, in the product: the per-trust approval toggle.</div>
            </div>
            <div class="panel" style="margin-top: 0.4em">
              <h3>Assessed against</h3>
              <p class="small" style="margin-bottom: 0">
                <span class="pill">Cyber Essentials</span><span class="pill">NHS DSPT / NDG</span
                ><span class="pill">NCSC Cyber Assessment Framework</span
                ><span class="pill">SATRE</span><span class="pill">Five Safes</span>
              </p>
            </div>
          </div>
        </div>
        <p class="small muted" style="margin-top: 0.25em; margin-bottom: 0; text-align: center">
          <strong>Patient data never leaves the hospital.</strong> Every layer between a researcher
          and a patient record enforces that boundary.
        </p>
      </div>
      <aside class="notes">
        (~1 min 30) (at 8:05) This is the slide your NHS and security audience came for, so give
        the compliance names their own beat — Cyber Essentials, DSPT/NDG, NCSC CAF, SATRE and the
        Five Safes are the words their IG teams use. Then the six guarantees. Guarantee 03 is the
        one worth pausing on: it is structural rather than a policy promise, which means a Trust
        does not have to trust our access controls at all, because there is no path in. If someone
        asks how the approval gate can be enforced when the hub is ours, that is the answer — the
        gate lives on the Trust's side. Also honest to mention: the authority to set approval is
        currently one platform-wide permission, which we found by auditing our own code. It is on
        the list.
      </aside>
    </section>

    <!-- 10 · The live run + job types. Third woven demo slide, and the payoff:
         the actual UK-Thailand deployment training in one job. Cambridge spent
         three slides here (deployment map, live run, recorded walkthroughs);
         this compresses to one, with the job-type families as the generalising
         claim and the two walkthrough QRs as the takeaway for anyone who wants
         the long version. -->
    <section>
      <div class="eyebrow">The platform · a real run, UK ⇄ Thailand</div>
      <h2>A federated run, live, per site</h2>
      <div class="slide-body">
        <div class="fig-split" style="--cols: 1.15fr 1fr; align-items: center">
          <div class="figure center">
            <img
              :src="asset('flip-live-training-bdms-kcl.png')"
              alt="FLIP admin UI showing a live federated training run, with the activity feed naming King's College London and Bangkok Dusit Medical Services as the selected trusts, and per-site training-loss curves"
              style="width: 100%; height: auto"
            />
            <div class="caption">
              KCL and Bangkok Dusit Medical Services training in the same job — Ark+ fine-tuning on
              five CXR pathology classes, 50 FedAvg rounds at ~27 KB per client per round.
            </div>
          </div>
          <div>
            <p class="small">
              <strong>Ark+</strong> (Ma et al., 2025) was pre-trained on 704,363 chest X-rays from
              six public datasets. We freeze the backbone and train only the
              <strong>5-lesion head — 6,885 parameters</strong>, which makes the round cheap.
            </p>
            <p class="small">
              A <strong>job type</strong> fixes the server-side contract: orchestration, aggregation
              and checkpoint staging. The researcher's <strong>app bundle</strong> fills the client
              side — a new app runs with <strong>no Central Hub or Trust Node code changes</strong>.
            </p>
            <div class="cols compact" style="--n: 2; gap: 0.5em; margin-top: 0.35em">
              <div class="panel center small">
                <h3>Classification</h3>
                <p class="small muted">CXR, Ark+ fine-tuning.</p>
                <span class="pill">FLARE</span><span class="pill">Flower</span>
              </div>
              <div class="panel center small">
                <h3>Segmentation</h3>
                <p class="small muted">3D spleen, MONAI bundle.</p>
                <span class="pill">FLARE</span><span class="pill">Flower</span>
              </div>
              <div class="panel center small">
                <h3>Evaluation</h3>
                <p class="small muted">Benchmarking, DeLong tests.</p>
                <span class="pill">FLARE</span><span class="pill">Flower</span>
              </div>
              <div class="panel center small">
                <h3>Synthesis</h3>
                <p class="small muted">Latent diffusion, federated.</p>
                <span class="pill">FLARE</span>
              </div>
            </div>
            <div style="display: flex; gap: 0.6em; margin-top: 0.4em; align-items: center">
              <img
                :src="reusableAsset('demo_video_xray.png')"
                alt="QR code linking to youtube.com/watch?v=BH97Yw_lmpw, the FLIP chest X-ray classification walkthrough"
                style="height: 70px; width: 70px; border-radius: 6px; flex: none"
              />
              <img
                :src="reusableAsset('demo_video_ct.png')"
                alt="QR code linking to youtube.com/watch?v=h7fHHkWEuEA, the FLIP CT spleen segmentation walkthrough"
                style="height: 70px; width: 70px; border-radius: 6px; flex: none"
              />
              <div class="caption" style="text-align: left; margin-top: 0">
                Full walkthroughs:<br>2D CXR (Ark+)<br>3D spleen (MONAI bundle)
              </div>
            </div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min 30) (at 9:20) The deployment is KCL and BDMS, a private Thai hospital
        group — cross-continental and cross-sector, which is a harder case than
        two NHS trusts. The generalising claim is the job-type grid: 2D
        classification and 3D organ segmentation share no application code, and
        neither needed a platform change. That is the reusability argument, and the
        MONAI bundle is the proof. Do not play the videos live — the QRs are there
        so the room can take them away.
      </aside>
    </section>

    <!-- 11 · NEW · A live clinical study. From the NECTEC deck's §14. The deck
         otherwise shows the platform and one proof-of-concept run but no
         standing clinical study, which is the difference between "this works"
         and "this is in use". The NECTEC slide carried UCLH and Guy's/King's as
         logos (not reusable), so the site list is written as text.
         UCLH is called out deliberately: it is TRE-based rather than running
         its own enclave, which is a genuine test of the deployment model, and
         it is the detail a platform audience will ask about. -->
    <section>
      <div class="eyebrow">In clinical use · a live study</div>
      <h2>FLIP-SIDE: segmenting organs for deep evaluation</h2>
      <div class="slide-body">
        <div class="cols" style="--n: 2; gap: 0.8em">
          <div class="panel flip">
            <h3>Primary objective</h3>
            <p class="small">
              Build an AI model that accurately segments <strong>internal organs from CT</strong>
              and extracts <strong>radiology biomarkers</strong> — trained federated across
              Guy's and St Thomas', King's College Hospital and UCLH.
            </p>
          </div>
          <div class="panel">
            <h3>Secondary objective</h3>
            <p class="small">
              Explore the role of radiology biomarkers <strong>combined with EHR data</strong>
              towards cancer treatment toxicity and acute illness outcomes.
            </p>
          </div>
        </div>
        <div class="cols compact" style="--n: 3; margin-top: 0.5em">
          <div class="panel center small">
            <h3>Three Trusts, one model</h3>
            <p class="small muted">Multi-site by default — no site holds the whole cohort.</p>
          </div>
          <div class="panel center small">
            <h3>Mixed deployment models</h3>
            <p class="small muted">UCLH joins as a <strong>TRE</strong> rather than running its own
              enclave — a useful test that the node is not tied to one hosting pattern.</p>
          </div>
          <div class="panel center small">
            <h3>Imaging + EHR together</h3>
            <p class="small muted">The multi-modal case neither the imaging estate nor the text
              estate can answer alone.</p>
          </div>
        </div>
        <p class="small muted" style="margin-top: 0.4em">
          Segmentation is a MONAI bundle job type; the biomarkers it produces become structured
          inputs to the next federated model.
        </p>
      </div>
      <aside class="notes">
        (~1 min) (at 11:00) This is the slide that shows FLIP is not a demo. Someone in this room
        will ask "who is actually using it" — FLIP-SIDE is the answer, and it is multi-site in
        production rather than a pilot. The UCLH point is worth a sentence because this audience
        cares about it: UCLH is TRE-based rather than a trust-run enclave, and FLIP's node ran
        there without a platform change, which is the strongest evidence that the deployment model
        generalises. The secondary objective is where the interesting science is — imaging
        biomarkers plus EHR to predict treatment toxicity.
      </aside>
    </section>

    <!-- 12 · NEW · Where the hub runs, and what is next. Merges the NECTEC
         deck's §6 (AWS Landing Zone Accelerator) with §16 (next steps). The
         AWS architecture diagram is the one image from that deck worth reusing:
         it is a real diagram, it shows the DGX A100s on-prem on the KCL side of
         the wire (which is the detail an NVIDIA audience looks for), and it
         makes the outbound-only claim concrete — the hub sits in a VPC and
         still has no inbound route to the Trust.
         Laid out as a full-width diagram over three short columns rather than
         a two-column split: the right-hand column carried two four-item lists
         and pushed the slide ~340px past the 720px box. -->
    <section>
      <div class="eyebrow">Where the hub runs · and what's next</div>
      <h2>Hub on AWS, nodes on site</h2>
      <div class="slide-body">
        <div class="figure center">
          <img
            :src="deckAsset('flip-aws-lza.png')"
            alt="AWS architecture diagram for FLIP: a networking account holding the VPC, route table, transit gateway and attachment; a FLIP production account with Central Hub DB on RDS Postgres, Central Hub APIs on ECS, a network load balancer, application load balancer and NAT gateway in public and private subnets; below, on-premise secure enclaves at KCL with a DGX A100 running ten 40GB GPUs, connected outbound to the hub"
            style="width: 100%; height: auto; max-height: 150px; object-fit: contain"
          />
          <div class="caption">
            The half that matters is the bottom: <strong>DGX A100</strong> inside the trust,
            only crossing link outbound.
          </div>
        </div>
        <div class="cols" style="--n: 4; gap: 0.55em; margin-top: 0.4em; align-items: start">
          <div class="panel flip">
            <h3>AWS-hosted hub</h3>
            <p class="small muted" style="margin-bottom: 0">
              CloudFront UI, users via <strong>Cognito</strong>. <strong>RDS</strong> and
              <strong>ECS</strong> underneath.
            </p>
          </div>
          <div class="panel">
            <h3>Landing Zone Accelerator</h3>
            <p class="small muted" style="margin-bottom: 0">
              Networking account for ingress, workload account for the hub. Sites are
              <strong>outbound only</strong>.
            </p>
          </div>
          <div class="panel">
            <h3>Next: multi-cloud nodes</h3>
            <p class="small muted" style="margin-bottom: 0">
              IaC for <strong>AWS, Azure, Snowflake</strong>; nodes in
              <strong>Leeds, Cambridge, Southampton</strong>.
            </p>
          </div>
          <div class="panel">
            <h3>Next: in-node labelling</h3>
            <p class="small muted" style="margin-bottom: 0">
              <strong>LLM auto-labelling</strong> inside the node — labels never leave either.
            </p>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min 15) (at 12:00) Two things. First, the hub is ordinary cloud infrastructure — the
        interesting part of that diagram is the bottom half, where the DGX sits inside the trust
        and the only arrow crossing the boundary is outbound. Second, next steps: the deployment
        shape is becoming a config, not an engineering project, which is what makes four more
        Trusts realistic in six months. If anyone wants to deploy a node, that is the conversation
        to have at the break — say it explicitly, this is a recruiting audience as much as a
        technical one.
      </aside>
    </section>

    <!-- 13 · Lessons learned. NEW for this deck, and promised by the session
         abstract: "real-world lessons learned in managing certificate
         provisioning, data heterogeneity, and regulatory hurdles". Cambridge
         had no equivalent slide. This is the slide a practitioner audience will
         actually use, so it gets plain, specific failures rather than
         platitudes. -->
    <section>
      <div class="eyebrow">What we got wrong first</div>
      <h2>Three lessons that cost us months</h2>
      <div class="slide-body">
        <div class="cols" style="--n: 3">
          <div class="panel">
            <h3>Certificates expire on a Sunday</h3>
            <p class="small">
              Startup kits are minted per site with a fixed lifetime, and a site that cannot rotate
              without a human is a site that silently drops out of the federation. Rotation had to
              become a first-class operation with its own alerting, not a provisioning afterthought.
            </p>
          </div>
          <div class="panel">
            <h3>Heterogeneity is not statistical</h3>
            <p class="small">
              The textbook worry is non-IID label distributions. The real one is that the same
              procedure is coded differently at each site, and a scanner at one trust writes a DICOM
              tag another leaves empty. <strong>OMOP + MI-CDM</strong> is where we spend the effort,
              because the harmonisation problem arrives before the training problem.
            </p>
          </div>
          <div class="panel flip">
            <h3>Governance is the long pole</h3>
            <p class="small">
              Technical onboarding of a new trust is days. Information governance approval is months,
              and it runs per site and per project. Designing the platform so each site's veto is
              enforced in code is what makes that conversation finishable at all.
            </p>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min 15) (at 10:50) The abstract promised these three, so deliver them
        plainly. The certificate one usually gets rueful laughter from anyone who
        has run a FLARE deployment for more than a year. The heterogeneity point is
        the one worth pressing: this audience's instinct is to reach for an
        algorithmic fix for non-IID data, and our experience is that the expensive
        problem is upstream of that and is a data-modelling problem.
      </aside>
    </section>

    <!-- 14 · The ecosystem, fused with the comparison and re-dotted for THIS
         event. Replaces Cambridge's slide 19 wholesale.
         Cambridge dotted platforms whose authors were at St John's. That set is
         irrelevant here, so every .here dot below was re-checked against the
         published FLARE Day 2026 programme (events.nvidia.com/flare-day-2026,
         read 14 Sep 2026) and maps to a talk in this event:
           NVIDIA FLARE  → the host; Cnudde, Chen, Xu, Roth, Belgiovine
           Mayo Clinic Platform → Cong (Charlie) Qin, 8:45 Sep 16
           Rhino          → Noy Maimon, 11:10 Sep 16 (shared-pool architecture)
           Apheris        → Nicolas Gautier, 9:05 Sep 16 (federated OpenFold3)
           Duality        → Omer Moran, 11:30 Sep 16
           CAIA           → Brian M. Bot, 8:05 Sep 16
           FLAIMME / NCI  → Umit Topaloglu, 12:10 Sep 16
           AIRR / BriCS   → Paul Wright, 10:30 Sep 16 — the talk immediately
                            after this one, hence the callout box.
         Platforms without a dot (CODA, Kaapana, MedPerf, FeatureCloud,
         vantage6, DataSHIELD, Swarm Learning, OpenFL, Substra, FATE, FedML,
         Fed-BioMed, FEDn) are from our comparison study and are NOT at this
         event — do not imply otherwise from the stage.
         The capability strip at the bottom is the short cut of
         PlatformComparisonDeck.vue: four columns instead of nine, chosen because
         they are the ones where the deployed-healthcare peer group is
         near-empty and where FLIP's position is therefore actually
         informative. -->
    <section>
      <div class="eyebrow">The wider ecosystem · and who is in this event</div>
      <h2>An ecosystem, not a race</h2>
      <div class="slide-body">
        <p class="small" style="margin: 0 0 0.3em">
          We audited 30 federated platforms for a MICCAI/DeCaF paper, scoring every cell from a local
          clone rather than from the paper. Cohort query and PACS are nearly empty across the peer
          group — that combination is where FLIP fits.
        </p>
        <div class="fig-split" style="--cols: 1.05fr 1fr; gap: 0.8em; align-items: start; margin-top: 0">
          <div>
            <table class="cap-table" style="font-size: 0.4em">
              <colgroup>
                <col style="width: 30%" />
                <col style="width: 26%" />
                <col span="4" style="width: 11%" />
              </colgroup>
              <thead>
                <tr>
                  <th>Platform</th>
                  <th>Data model</th>
                  <th>Cohort</th>
                  <th>PACS</th>
                  <th>Appr.</th>
                  <th>Prov.</th>
                </tr>
              </thead>
              <tbody>
                <tr class="is-ours">
                  <td><strong>FLIP (ours)</strong></td>
                  <td>OMOP + MI-CDM</td>
                  <td><span class="cov yes">●</span></td>
                  <td><span class="cov yes">●</span></td>
                  <td><span class="cov partial">◐</span></td>
                  <td><span class="cov yes">●</span></td>
                </tr>
                <tr>
                  <td>CODA</td>
                  <td>FHIR + DICOM</td>
                  <td><span class="cov yes">●</span></td>
                  <td><span class="cov partial">◐</span></td>
                  <td><span class="cov no">✕</span></td>
                  <td><span class="cov partial">◐</span></td>
                </tr>
                <tr>
                  <td>JIP / Kaapana</td>
                  <td>DICOM</td>
                  <td><span class="cov partial">◐</span></td>
                  <td><span class="cov yes">●</span></td>
                  <td><span class="cov yes">●</span></td>
                  <td><span class="cov partial">◐</span></td>
                </tr>
                <tr>
                  <td>MedPerf</td>
                  <td>Benchmark-defined</td>
                  <td><span class="cov partial">◐</span></td>
                  <td><span class="cov no">✕</span></td>
                  <td><span class="cov yes">●</span></td>
                  <td><span class="cov partial">◐</span></td>
                </tr>
                <tr>
                  <td>FeTS</td>
                  <td>NIfTI / BraTS</td>
                  <td><span class="cov no">✕</span></td>
                  <td><span class="cov no">✕</span></td>
                  <td><span class="cov no">✕</span></td>
                  <td><span class="cov no">✕</span></td>
                </tr>
              </tbody>
            </table>
            <p class="small muted" style="margin-top: 0.25em; font-size: 0.5em; margin-bottom: 0">
              <span class="cov yes">●</span>ships &nbsp;<span class="cov partial">◐</span>partial
              &nbsp;<span class="cov no">✕</span>no. 5 of 30 audited shown. Kaapana beats us on
              approval.
            </p>
            <div class="flip-note" style="margin-top: 0.4em">
              <span class="flip-tag">Next talk</span> Paul Wright follows with federated learning on
              the UK's AIRR and Isambard-AI. FLIP runs where the data is; AIRR is where the compute
              is — two ends of the same UK problem.
            </div>
          </div>
          <div>
            <div class="panel flip">
              <h3>Also solving this, and in this event</h3>
              <div style="margin-top: 0.2em">
                <span class="pill ours here">FLIP</span><span class="pill here">Mayo Clinic Platform</span
                ><span class="pill here">Cancer AI Alliance</span><span class="pill here">FLAIMME · NCI</span
                ><span class="pill here">Apheris</span><span class="pill here">Rhino</span
                ><span class="pill here">Duality</span><span class="pill here">AIRR · BriCS</span>
              </div>
              <p class="small muted" style="margin: 0.4em 0 0; font-size: 0.62em">
                <span class="here-dot"></span>presenting at FLARE Day 2026. Undotted platforms in the
                table are from the audit, not the agenda.
              </p>
            </div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~1 min 30) (at 12:05) Two jobs on this slide. First, the honest version of
        the comparison: the cohort-query and PACS columns are nearly empty across
        the deployed peer group, and that combination is where FLIP fits — a
        statement about orchestration, not about being better. Kaapana beats us on
        per-site approval and its runtime policy enforcement is genuinely stronger
        than ours (Traefik forwardAuth in front of an Open Policy Agent decision,
        blocking the request path rather than the rendered page). Say that out loud;
        it is why the table reads as evidence rather than marketing. Second, the
        room: everything with a dot is presenting here today. Explicitly hand over
        to Paul Wright — his talk is next and it is the complementary half of ours.
        Careful not to imply the undotted platforms are here; they are from the
        audit, not the agenda.
      </aside>
    </section>

    <!-- 15 · The team. Kept from Cambridge, trimmed to one row of institutions.
         This deck is billed as a three-institution open-source project, so the
         GSTT column is doing real work in the credibility argument. -->
    <section>
      <div class="eyebrow">Building FLIP together</div>
      <h2>FLIP team</h2>
      <div class="slide-body">
        <div class="people" style="margin-top: 0">
          <div class="photo-chip">
            <img
              :src="sharedAsset('shared/seb_kcl.png')"
              alt="Professor Sébastien Ourselin, Head of School, School of Biomedical Engineering & Imaging Sciences"
              style="max-height: 150px"
            />
          </div>
        </div>
        <div class="fig-split" style="--cols: 4fr 3fr; margin-top: 0.3em; gap: 1.2em; align-items: start">
          <div>
            <div class="eyebrow" style="text-align: center; color: #eb2f2d; font-weight: 600; font-size: large;">King's College London</div>
            <div class="people">
              <div class="photo-chip">
                <img :src="sharedAsset('shared/jorge.png')" alt="Dr M. Jorge Cardoso, Group Lead & Reader" style="height: 155px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/rafa_kcl.png')" alt="Rafael Garcia-Dias, Senior AI Engineer on Foundational Models for Healthcare" style="height: 155px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/alex_kcl.png')" alt="Alexandre Triay Bagur, Senior AI Engineer" style="height: 155px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/virginia_kcl.png')" alt="Virginia Fernandez, Research Associate" style="height: 155px" />
              </div>
            </div>
          </div>
          <div>
            <div class="eyebrow" style="text-align: center; color: #005EB8; font-weight: 600; font-size: large;">Guy's and St Thomas' Trust</div>
            <div class="people">
              <div class="photo-chip">
                <img :src="sharedAsset('shared/joe_gstt.png')" alt="Joe Zhang, Head of Data Science" style="height: 155px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/lawrence_gstt.png')" alt="Lawrence Adams, Lead Analytics Engineer" style="height: 155px" />
              </div>
              <div class="photo-chip">
                <img :src="sharedAsset('shared/martin_gstt.png')" alt="Martin Chapman, Lead NLP Engineer" style="height: 155px" />
              </div>
            </div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~20s) (at 13:35) Quick. The point is that engineering, clinical and
        governance people sit in one group, and that is itself part of why the
        deployment happened.
      </aside>
    </section>

    <!-- 16 · Close. Stays up through Q&A so the QR codes remain on screen. -->
    <section class="title-slide center">
      <div class="eyebrow">Federated Learning &amp; Interoperability Platform</div>
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
            <p class="figure-placeholder__desc" style="margin-bottom: 0.35em">Live demo</p>
            <div class="figure-placeholder compact" style="aspect-ratio: 1 / 1; min-height: 0; padding: 0.55em">
              <img
                :src="reusableAsset('ark_demo.png')"
                alt="https://app.flip.aicentre.co.uk/ark_demo/"
                title="A read-only snapshot of the Ark+ federated experiments in the real FLIP interface"
                style="display: block; width: 100%; height: 100%; object-fit: cover; border-radius: 6px"
              />
            </div>
          </div>
        </div>
        <p class="venue-note" style="margin-top: 0.6em">
          The third QR is a read-only snapshot of the UK ⇄ Thailand Ark+ runs, in the real interface.
        </p>
      </div>
      <aside class="notes">
        (~20s, then Q&amp;A to 10:30) (at 13:55) Leave this up. Three QRs: contact,
        the repository, and the live read-only demo — that last one is the one this
        audience will actually scan, so say what it is.
      </aside>
    </section>
  </RevealDeck>
</template>
