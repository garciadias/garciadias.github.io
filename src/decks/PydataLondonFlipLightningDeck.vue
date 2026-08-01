<script setup>
import DeckBrand from '@/components/DeckBrand.vue';
import RevealDeck from '@/components/RevealDeck.vue';

// A 10-slide lightning talk for a PyData London audience. Narrative: imaging
// is the largest untapped modality in UK health data (locked in PACS + free
// text reports) → two things block it at national scale, an infrastructure
// problem (bespoke federated nodes) and a trust problem (privacy isn't
// automatic in FL) → the result is one-off pilots, not standing
// infrastructure → FLIP fixes both, reusable + privacy-preserving + governed.
// Gif spots are left as dashed placeholders (see .figure-placeholder in
// RevealDeck.vue); drop a real file at the named path under
// public/presentations/pydata-london-flip-lightning-2026/ and swap the
// placeholder <div> for an <img>/<video> once you've picked the gifs.
const asset = (name) => `${import.meta.env.BASE_URL}presentations/pydata-london-flip-lightning-2026/${name}`
const sharedAsset = (name) => `${import.meta.env.BASE_URL}presentations/${name}`
</script>

<template>
  <RevealDeck :options="{ center: true }">
    <template #chrome>
      <DeckBrand />
    </template>

    <!-- 0 · Title -->
    <section class="title-slide center">
      <div class="eyebrow">PyData London ⚡ Lightning talk ⚡ FLIP: Federated Learning and Interoperability Platform</div>
      <h1>AI on the NHS: can federated learning preserve patient privacy?</h1>
      <p class="subtitle">
        Yes, but only if privacy is architecture, not a promise. Here's how FLIP unlocks the
        millions of NHS scans, without exposing patient data.
      </p>
      <p class="venue-note">Rafael Garcia-Dias · AI Centre for Value-Based Healthcare · King's College London</p>
      <aside class="notes">
        (~20s) Hook fast: imaging is the biggest untapped modality in UK health data, and the reason
        isn't the maths — it's engineering and trust. That's the whole talk in one sentence.
      </aside>
    </section>

    <!-- 1 · The opportunity -->
    <section>
      <div class="eyebrow">The opportunity</div>
      <h2>Every major centre is sitting on millions of studies</h2>
      <div class="fig-split" style="--cols: 1.35fr 0.9fr; margin-top: 0.3em; align-items: start">
        <div>
          <p class="small">
            Every major NHS trust has a PACS archive going back decades. <br/>
            It holds X-rays, CTs, and MRIs, and every single one already has an expert interpretation attached to it.
          </p>
          <ul class="checklist small">
            <li>Millions of studies per trust, accumulating every year</li>
            <li>Almost none of it has ever trained a model at national scale</li>
          </ul>
        </div>
        <div class="figure" style="aspect-ratio: 4 / 3">
          <img
            :src="asset('nhs_medical_images.png')"
            alt="Line chart of annual NHS medical imaging procedures in England, rising from 39.94M in 2014/15 to 47.16M in 2023/24, with a dip to 34.92M during the 2020/21 COVID lockdown"
          />
        </div>
      </div>
      <aside class="notes">
        (~25s) Open with scale, not the problem — this isn't a small-data problem, it's a huge,
        already-labelled dataset that nobody can train on yet.
      </aside>
    </section>

    <!-- 2 · Problem one: infrastructure -->
    <section>
      <div class="eyebrow">Problem one · infrastructure</div>
      <h2>FLIP it! Sending the model to the data</h2>
        <p class="small">
          Imaging data is too large and too sensitive to centralise. Law, policy, governance, and ethics block the transfer of raw data, so the model has to go to the data, not the other way around.
        </p>
      <div class="fig-split" style="--cols: 1.4fr 0.8fr; margin-top: 0.3em">
        <div>
        <p class="small">
          But making it part of a permanent, reusable platform is another story.
        </p>
          <div class="cols" style="--n: 2; margin-top: 0.6em">
            <div class="panel">
              <h3>Requirements</h3>
              <ul class="checklist small">
                <li>Data never leaves the site</li>
                <li>Training happens locally</li>
                <li>No PII in model assets</li>
              </ul>
            </div>
          </div>
        </div>
      </div>
      <aside class="notes">
        (~25s) This is the infrastructure half of the problem — federating isn't hard conceptually,
        standing up the plumbing at every site, repeatedly, is.
      </aside>
    </section>

    <!-- 3 · Problem one: infrastructure (cont.) -->
    <section data-transition="fade">
      <div class="eyebrow">Problem one · infrastructure</div>
      <h2>FLIP it! Sending the model to the data</h2>
        <p class="small">
          Imaging data is too large and too sensitive to centralise. Law, policy, governance, and ethics block the transfer of raw data, so the model has to go to the data, not the other way around.
        </p>
      <div class="fig-split" style="--cols: 1.4fr 0.8fr; margin-top: 0.3em">
        <div>
        <p class="small">
          But making it part of a permanent, reusable platform is another story.
        </p>
          <div class="cols" style="--n: 2; margin-top: 0.6em">
            <div class="panel">
              <h3>Requirements</h3>
              <ul class="checklist small">
                <li>Data never leaves the site</li>
                <li>Training happens locally</li>
                <li>No PII in model assets</li>
              </ul>
            </div>
            <div class="panel flip">
              <h3>Reality</h3>
              <ul class="crosslist small">
                <li>No standard node setup</li>
                <li>FL stacks are complex</li>
                <li>Short pilots, no budget</li>
              </ul>
            </div>
          </div>
        </div>
        <div class="figure" style="aspect-ratio: 5 / 3">
          <img
            :src="asset('no_time.gif')"
            alt="GIF of the 'ain't nobody got time for that' meme"
          />
        </div>
      </div>
      <aside class="notes">
        (~25s) This is the infrastructure half of the problem — federating isn't hard conceptually,
        standing up the plumbing at every site, repeatedly, is.
      </aside>
    </section>

    <!-- 4 · The result -->
    <section>
      <div class="eyebrow">The result</div>
      <h2>Going to production in real hospitals is challenging</h2>
      <div class="fig-split" style="--cols: 1.1fr 0.65fr; margin-top: 0.3em; align-items: start">
        <div>
          <ul class="checklist small">
            <li>One brilliant pilot presented once at one conference</li>
          </ul>
          <ul class="crosslist small">
            <li>"privacy-preserving" was never independently guaranteed</li>
            <li>No institution approves it for live infrastructure</li>
            <li>No shared infrastructure for the next hospital to plug into, so it starts from zero, again</li>
          </ul>
        </div>
        <div class="figure" style="width: fit-content; max-width: 100%; margin: 0 auto; text-align: center">
          <a href="https://imgflip.com/i/axwkg1" target="_blank" rel="noopener">
            <img
              src="https://i.imgflip.com/axwkg1.jpg"
              alt="Meme illustrating that federated imaging AI often ends up as one-off pilots without shared infrastructure or guaranteed privacy"
              style="max-height: 340px; max-width: 100%"
            />
          </a>
        </div>
      </div>
      <aside class="notes">
        (~20s) Land this as the punchline of the problem section — pause here. This is the "so what"
        the rest of the talk answers.
      </aside>
    </section>

    <!-- 5 · Problem two: trust -->
    <section>
      <div class="eyebrow">Problem two · trust</div>
      <h2>Your phone is doing it and you should be concerned</h2>
          <p class="small">
            Google Gboard is the clearest, best-documented case: researchers showed that the federated
            system used to train next-word prediction is not private.
          </p>
      <div class="fig-split" style="--cols: 1.15fr 1fr; margin-top: 0.3em">
        <div>
          <ul class="crosslist small">
            <li>Model updates could be attacked to reconstruct the actual sentences a user typed</li>
            <li>Common countermeasures like adding noise only helped if they damaged model utility</li>
            <li>That directly contradicts the privacy promise of FL in a production app with billions of downloads</li>
          </ul>
          <p class="small muted" style="margin-top: 0.5em">
            Source: <a href="https://arxiv.org/abs/2210.16947" target="_blank" rel="noopener">Suliman and Leith, 2023</a>
          </p>
        </div>
        <div class="figure" style="width: fit-content; max-width: 100%; margin: 0 auto; text-align: center">
          <img
            src="https://media2.giphy.com/media/v1.Y2lkPTc5MGI3NjExcW5oYnN3aW8yZzZocHl1a2xlazA1ejkzcGlja3FkZnlyNG5pM3J5biZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/WAnUdXKHObQUU/giphy.gif"
            alt="Animated GIF emphasizing privacy fine print and unexpected leakage"
            style="display: block; width: auto; height: auto; max-width: 100%; max-height: 100%; margin: 0 auto"
          />
        </div>
      </div>
      <aside class="notes">
        (~25s) Use Gboard as the evidence slide: this is not a theoretical concern, and the privacy
        failure survives standard fixes unless they make the model less useful.
      </aside>
    </section>

    <!-- 6 · Enter FLIP -->
    <section>
      <div class="eyebrow">Enter FLIP</div>
      <h2>FLIP turns a one-off PoC into a reusable platform</h2>
      <p class="small">
        FLIP is a multi-application, open-source federated learning platform. It packages the work
        needed for medical imaging studies into a repeatable layer, so each new project can reuse the
        same infrastructure instead of rebuilding it.
      </p>
      <div class="cols" style="--n: 3; margin-top: 0.7em">
        <div class="panel">
          <h3>Reusable</h3>
          <p class="small">The same platform supports different FL jobs, from cohort-driven training to federated evaluation.</p>
        </div>
        <div class="panel">
          <h3>Lifecycle</h3>
          <p class="small">It covers the full project flow: cohort discovery, staging, local execution, secure aggregation, and provenance.</p>
        </div>
        <div class="panel flip">
          <h3>Built for imaging</h3>
          <p class="small">It already powers medical-imaging use cases, but it is not limited to one modality or one study.</p>
        </div>
      </div>
      <aside class="notes">
        (~25s) Frame FLIP as the platform layer: reusable across applications, not tied to one study,
        and responsible for the mechanics that usually get rebuilt every time.
      </aside>
    </section>

    <!-- 7 · How it works -->
    <section>
      <div class="eyebrow">How it works</div>
      <h2>FLIP packages the full federated workflow</h2>
      <div class="fig-split" style="--cols: 1fr 2.18fr; margin-top: 0.5em; align-items: center">
        <ol class="contribs small">
          <li><strong>Define once.</strong> Write the cohort query, upload the app bundle, and pick the job type.</li>
          <li><strong>Approve locally.</strong> Each institution reviews the project through its own governance process.</li>
          <li><strong>Run locally.</strong> The Trust Node runs the job and returns only weights, checkpoints, and results.</li>
        </ol>
        <div class="figure center">
          <img
            :src="sharedAsset('flip-maturity-pitch-2026/flip-architecture-uk-thailand.png')"
            alt="Symmetric FLIP architecture: two mirrored trust nodes, each taking PACS, EHR, and other clinical systems into a secure enclave holding a structured OMOP database and XNAT, feeding an AI training client, connected through a FLIP node to the platform"
            style="max-height: 290px; max-width: 100%"
          />
        </div>
      </div>
      <aside class="notes">
        (~25s) Walk left to right: define the cohort once, get site-level approval, then run locally.
        The point is the repeatable workflow, not the diagram labels.
      </aside>
    </section>

    <!-- 8 · Governance, revisited -->
    <section>
      <div class="eyebrow">Back to problem two</div>
      <h2>Privacy here is architecture, not a promise</h2>
      <ul class="checklist small">
        <li>Only model weights, checkpoints, and aggregate stats leave the hospital</li>
        <li>Images and row-level patient data never leave the Trust Node</li>
        <li>The NHS National Data Opt-Out is honoured automatically, by design</li>
        <li>Every site can veto any project independently through its own governance process</li>
        <li>Every round is logged and auditable, so provenance is built into the platform</li>
      </ul>
      <aside class="notes">
        (~25s) This is the direct answer to the trust problem: the platform enforces the boundary,
        while governance decides whether the project runs at all.
      </aside>
    </section>

    <!-- 9 · Thank you -->
    <section class="title-slide center">
      <div class="eyebrow">Federated Learning and Interoperability Platform</div>
      <h1>Contribute: Help build FLIP</h1>
      <p class="byline"><a href="https://github.com/londonaicentre/FLIP">github.com/londonaicentre/FLIP</a> · Apache 2.0</p>
      <div class="cols" style="--n: 3; gap: 0.7em; max-width: 22em; margin: 0.8em auto 0; align-items: stretch">
        <div>
          <p class="figure-placeholder__desc" style="margin-bottom: 0.35em">Contact me</p>
          <div class="figure-placeholder compact" style="aspect-ratio: 1 / 1; min-height: 0; padding: 0.55em">
            <img
              :src="asset('me.png')"
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
              :src="asset('FLIP_REPO.png')"
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
              :src="asset('ai_centre.png')"
              alt="https://www.aicentre.co.uk/about/meet-the-team"
              title="AI Centre team"
              style="display: block; width: 100%; height: 100%; object-fit: cover; border-radius: 6px"
            />
          </div>
        </div>
      </div>
      <aside class="notes">
        (~10s) Close on the repo: invite people to try it, file issues, and contribute code or feedback.
      </aside>
    </section>
  </RevealDeck>
</template>
