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
    id: 'iaa-so-chemical-tagging-2026',
    title: 'Chemical tagging: finding lost star clusters',
    subtitle: 'From K-means to EVoC — clustering chemistry to reconstruct dissolved clusters',
    date: '2026',
    venue: 'IAA-SO School on AI/ML in Astronomy 2026 · Unsupervised Learning',
    cover: `${import.meta.env.BASE_URL}presentations/iaa-so-chemical-tagging-2026/benchmark_grid.png`,
    description:
      'A ~90-minute lecture for the IAA-SO School on AI/ML in Astronomy (unsupervised-learning pillar). ' +
      'It walks the clustering toolbox one algorithm at a time, each chosen to fix the failure of the last — ' +
      'K-means, KNN, DBSCAN, HDBSCAN*, PLSCAN (Bot, McInnes & Aerts 2025), t-SNE, UMAP, and finally EVoC, ' +
      'which fuses all three ideas — embedding, density clustering, and persistence-based scale selection — ' +
      'into one pass. The running example is chemical tagging: reconstructing dissolved star clusters from ' +
      '16-D APOGEE DR17 + Gaia EDR3 abundances (Kos et al. 2017), with every method scored against kinematic ' +
      'ground truth. Through-line: there are no bad models, only models applied outside their assumptions — ' +
      'so know your data and validate against reality.',
    tags: ['Astronomy', 'Chemical Tagging', 'Clustering', 'KNN', 't-SNE', 'UMAP', 'HDBSCAN', 'PLSCAN', 'EVoC', 'Unsupervised Learning'],
    deck: () => import('@/decks/IaaSoChemicalTaggingDeck.vue')
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
      'layer between the hospital and the framework, solved once, open source and Apache 2.0. Covers the ' +
      'services it decomposes into, the single outbound connection it asks a hospital network for, and ' +
      'where the guardrails execute — then a demo of a real UK ⇄ Thailand deployment: cohort query, ' +
      'per-site approval, and a live federated training run.',
    tags: ['Federated Learning', 'FLIP', 'NHS', 'Open Source', 'Platform', 'Cambridge'],
    deck: () => import('@/decks/FlipReusableFlPlatformDeck.vue')
  }
]

// Decks shown on the public gallery (everything not flagged `unlisted`).
export const listedPresentations = presentations.filter((p) => !p.unlisted)

export function getPresentation(id) {
  return presentations.find((p) => p.id === id)
}
