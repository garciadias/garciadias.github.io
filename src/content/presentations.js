// Registry of slide decks. Each entry feeds the /presentations gallery and the
// fullscreen /presentations/:id deck route. `deck` is a lazy import of the Vue
// component holding the reveal.js <section> slides.
//
// To add a new presentation: create src/decks/<Name>Deck.vue, then add an entry
// here. Nothing else needs editing.
export const presentations = [
  {
    id: 'fla3-governance-federated-learning',
    title: 'Governance as a First-Class Control',
    subtitle: 'FLA³ (Federated Learning with AAA) — read throughout against FLIP',
    date: '2026',
    venue: 'Paper walkthrough · arXiv 2603.10063',
    description:
      'A technical walkthrough of FLA³ — a federated learning platform that enforces Authentication, Authorisation and Accounting as a runtime control plane (XACML PDP, per-round, fail-closed, cryptographically audited). The argument: in regulated healthcare the breach is unauthorised computation, not data movement — and governance can be enforced at zero accuracy cost. FLIP is woven in throughout as a peer system.',
    tags: ['Federated Learning', 'Governance / AAA', 'XACML', 'Healthcare AI', 'FLIP'],
    deck: () => import('@/decks/Fla3Deck.vue')
  }
]

export function getPresentation(id) {
  return presentations.find((p) => p.id === id)
}
