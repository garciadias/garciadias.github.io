// Registry of slide decks. Each entry feeds the /presentations gallery and the
// fullscreen /presentations/:id deck route. `deck` is a lazy import of the Vue
// component holding the reveal.js <section> slides.
//
// To add a new presentation: create src/decks/<Name>Deck.vue, then add an entry
// here. Nothing else needs editing.
//
// Set `unlisted: true` to publish a deck without showing it on the public
// /presentations gallery. Unlisted decks still work at their direct
// /presentations/:id URL and appear on the secret /presentations/list index.
export const presentations = [
  {
    id: 'fla3-governance-federated-learning',
    title: 'What can FLIP learn from FLA³?',
    subtitle: 'A FLIP-eyed read of FLA³ — Federated Learning with Authentication · Authorisation · Accounting',
    date: '2 Jun 2026',
    venue: 'AMIGO team meeting',
    cover: `${import.meta.env.BASE_URL}presentations/fla3/x4.png`,
    description:
      'A walkthrough of the FLA³ paper (arXiv 2603.10063) through one question: what can our FLIP platform take from it? FLA³ enforces Authentication, Authorisation & Accounting as a runtime control plane — per-round XACML policy decisions, fail-closed, with cryptographically signed audit. Its thesis: in regulated healthcare the breach is unauthorised computation, not data movement — and governance costs nothing in accuracy (federation matches centralised on the INTERVAL iron-deficiency task, and lifts the weakest sites most). Every slide carries an honest FLIP verdict: where we already match, and the upgrades worth weighing.',
    tags: ['Federated Learning', 'Governance', 'Healthcare AI', 'FLIP', 'FLA³'],
    deck: () => import('@/decks/Fla3Deck.vue')
  },
  {
    id: 'miccai_decaf_2026_draft',
    title: 'Federated chest X-ray learning, UK ⇄ Thailand',
    subtitle: 'Cross-continental federated fine-tuning of a CXR foundation model with FLIP — MICCAI 2026 / DECAF working draft',
    date: '2026',
    venue: 'MICCAI 2026 · DECAF (draft)',
    description:
      'A working-draft walkthrough of a UK–Thailand proof of concept: deploying FLIP — our Federated Learning & Interoperability Platform, validated at NHS scale — in a setting it was not designed for, a live cross-continental, cross-sector federation linking UK academia with a private Thai hospital group. Not a new algorithm: a reproducible federated fine-tuning workflow for a pretrained chest X-ray foundation model, a four-arm comparison (zero-shot · local-only · federated · centralised) over locally generated synthetic data, and an evaluation that pairs predictive metrics with system- and deployment-level measurements. Unlisted draft for internal review.',
    tags: ['Federated Learning', 'Chest X-ray', 'Foundation Models', 'FLIP', 'MICCAI 2026', 'Draft'],
    unlisted: true,
    deck: () => import('@/decks/MiccaiDecaf2026Deck.vue')
  },
  {
    id: 'pydata-london-flip-lightning-2026',
    title: 'AI on the NHS: can federated learning preserve patient privacy?',
    subtitle: 'A 5-minute case for sending the model to the data, not the other way round',
    date: '2026',
    venue: 'PyData London · Lightning talk',
    cover: `${import.meta.env.BASE_URL}presentations/flip-maturity-pitch-2026/flip-architecture-symmetric.png`,
    description:
      'Federated learning inverts the logic of big-data lakes for training AI models: instead of bringing ' +
      'data to models, we send the models to be trained where the data is. A fast, informal cut of FLIP — ' +
      'the platform the AI Centre for Value-Based Healthcare runs in production across the NHS and a ' +
      'Thai hospital group — built for a PyData lightning-talk slot: the problem, the one idea, the ' +
      'architecture in one diagram, real production numbers, and the open-source Python stack underneath it.',
    tags: ['Federated Learning', 'Privacy', 'NHS', 'PyData', 'Lightning Talk', 'FLIP'],
    deck: () => import('@/decks/PydataLondonFlipLightningDeck.vue')
  },
  {
    id: 'flip-maturity-pitch-2026',
    title: 'FLIP is production infrastructure now',
    subtitle: 'A year of engineering turned federated learning from a prototype into live clinical infrastructure — running today across the UK and Thailand.',
    date: 'Jul 2026',
    venue: 'AMIGO team meeting',
    cover: `${import.meta.env.BASE_URL}presentations/flip-maturity-pitch-2026/flip-architecture-symmetric.png`,
    description:
      'A platform status and roll-out briefing: five slides on the last year of FLIP engineering — from a ' +
      'single-machine prototype to live NHS and Thai trust nodes — then the shared case for going live, backed ' +
      'by dated, screenshotted evidence rather than roadmap promises.',
    tags: ['Federated Learning', 'FLIP', 'NHS', 'Platform', 'Production'],
    unlisted: true,
    deck: () => import('@/decks/FlipMaturityPitch2026Deck.vue')
  },
  {
    id: 'flip-reusable-fl-platform-2026',
    title: 'FLIP: a reusable, general-purpose open-source FL platform',
    subtitle: 'From Zero to Hero: Federated AI in Healthcare Systems · St John\'s College, Cambridge',
    date: '2026',
    venue: 'Federated AI in Healthcare Systems · Cambridge',
    cover: `${import.meta.env.BASE_URL}presentations/flip-reusable-fl-platform/intro-globe-poster.jpg`,
    description:
      'A 25-minute workshop talk on what it takes to turn federated learning from a research project ' +
      'into standing infrastructure. The argument: the FL frameworks are ready — Flower and NVIDIA FLARE ' +
      'made the federated core something you can build on — but every collaboration still rebuilds ' +
      'cohort definition, harmonisation, per-site approval, scheduling and audit by hand. FLIP is that ' +
      'layer between the hospital and the framework, solved once, open source and ' +
      'Apache 2.0. Centred on a demo of a real UK ⇄ Thailand deployment — cohort query, governance, live ' +
      'training run — followed by what the deployment measured, what actually broke, and the lessons ' +
      'that only show up in production.',
    tags: ['Federated Learning', 'FLIP', 'NHS', 'Open Source', 'Platform', 'Cambridge'],
    deck: () => import('@/decks/FlipReusableFlPlatformDeck.vue')
  },
  {
    id: 'flip-dev-report-2026-jun-jul',
    title: 'v0.3.0 shipped — one codebase, Kubernetes trusts, foundation-model FL',
    subtitle: 'FLIP development report · 1 Jun – 14 Jul 2026',
    date: '14 Jul 2026',
    venue: 'FLIP development report',
    cover: `${import.meta.env.BASE_URL}presentations/flip-dev-report-2026-jun-jul/slide-4k.png`,
    description:
      'Single-slide development report for 1 Jun – 14 Jul 2026: v0.3.0 shipped with one merged FL codebase ' +
      '(NVFLARE and Flower in-tree), Kubernetes trust deployment via Helm, and federated foundation-model ' +
      'fine-tuning of the Ark+ chest X-ray model — 105 PRs merged, 247k lines changed, 5 epics closed.',
    tags: ['FLIP', 'v0.3.0', 'Kubernetes', 'Foundation Models', 'Report'],
    unlisted: true,
    deck: () => import('@/decks/FlipReportJunJul2026Deck.vue')
  },
  {
    id: 'miccai-decaf-2026-results',
    title: 'Making Cross-Continental Federated Learning Repeatable with FLIP',
    subtitle: 'A multi-application study — federated fine-tuning & evaluation of a chest X-ray foundation model, UK ⇄ Thailand.',
    date: '28 Jul 2026',
    venue: 'AMIGO team meeting',
    cover: `${import.meta.env.BASE_URL}presentations/miccai-decaf-2026-results/flip-architecture-diagram.png`,
    description:
      'Results talk for the MICCAI 2026 / DECAF workshop submission (anonymised, under review): a ' +
      'multi-application study testing whether FLIP — reusable federated-learning infrastructure — holds up ' +
      'outside the NHS, exercising two distinct FL applications on real UK ⇄ Thailand sites.',
    tags: ['Federated Learning', 'Chest X-ray', 'Foundation Models', 'FLIP', 'MICCAI 2026', 'Results'],
    unlisted: true,
    deck: () => import('@/decks/MiccaiDecaf2026ResultsDeck.vue')
  },
  // The Boehringer interview deck lives in .data/ (gitignored) so it never
  // publishes to the public site. Its registry entry is parked here rather than
  // deleted: re-enable it only alongside a deck file that is safe to publish.
  //   To view it locally without committing it:
  //     ln -s ../../.data/BoehringerSeniorMLEngineerDeck.vue src/decks/
  //     echo 'src/decks/BoehringerSeniorMLEngineerDeck.vue' >> .gitignore
  //   then uncomment the block below.
  // {
  //   id: 'boehringer-senior-ml-engineer-2026',
  //   title: 'From Research Prototype to Production: Multimodal Foundation Models for Drug Discovery',
  //   subtitle: 'Senior ML Engineer · AI Accelerator · Boehringer Ingelheim',
  //   date: 'Aug 2026',
  //   venue: 'Interview presentation',
  //   description:
  //     'Interview presentation: taking a dual-encoder EHR + proteomics transformer from a research notebook ' +
  //     'to a production-grade, scalable, fine-tunable foundation model for drug discovery — without breaking ' +
  //     'the biology.',
  //   tags: ['Foundation Models', 'Drug Discovery', 'Multimodal', 'ML Engineering'],
  //   unlisted: true,
  //   deck: () => import('@/decks/BoehringerSeniorMLEngineerDeck.vue')
  // },
  // The Boehringer MLOps interview deck lives in .data/ (gitignored) so it never
  // publishes to the public site. Its registry entry is parked here rather than
  // deleted: re-enable it only alongside a deck file that is safe to publish.
  //   To view it locally without committing it:
  //     ln -s ../../.data/boehringer-mlops-engineer-2026/BoehringerSeniorMLOpsEngineerDeck.vue src/decks/
  //     echo 'src/decks/BoehringerSeniorMLOpsEngineerDeck.vue' >> .gitignore
  //   then uncomment the block below.
  {
    id: 'boehringer-senior-mlops-engineer-2026',
    title: 'Operationalising Foundation Models for Drug Discovery',
    subtitle: 'Senior MLOps Engineer · AI Enablement · Boehringer Ingelheim',
    date: '2026',
    venue: 'Interview presentation',
    description:
      'Interview presentation: an end-to-end MLOps operating model for training, serving and governing ' +
      'multimodal foundation models and their downstream fine-tuned models — lineage, reproducibility, ' +
      'compute cost, monitoring, handover and lifecycle.',
    tags: ['MLOps', 'Foundation Models', 'Drug Discovery', 'Lineage'],
    unlisted: true,
    deck: () => import('@/decks/BoehringerSeniorMLOpsEngineerDeck.vue')
  }
]

// Decks shown on the public gallery (everything not flagged `unlisted`).
export const listedPresentations = presentations.filter((p) => !p.unlisted)

export function getPresentation(id) {
  return presentations.find((p) => p.id === id)
}
