export const profile = {
  name: 'Dr Rafael Garcia-Dias',
  shortName: 'Rafael Garcia-Dias',
  headline: 'Senior AI Engineer — Federated Learning, Foundation Models & Healthcare AI at NHS Scale',
  location: 'London, UK',
  rightToWork: 'Right to work in the UK — no sponsorship required.',
  cvPath: '/static/img/cv/curriculum-vitae.pdf',
  socials: [
    { label: 'GitHub', href: 'https://github.com/garciadias', icon: 'github' },
    { label: 'LinkedIn', href: 'https://www.linkedin.com/in/garcia-dias/', icon: 'linkedin' },
    { label: 'Google Scholar', href: 'https://scholar.google.co.uk/citations?hl=en&user=MlwZerQAAAAJ', icon: 'scholar' }
  ],
  intro:
    'Senior AI Engineer with a PhD in Machine Learning and 10+ years of production Python experience. Currently building FLIP — the open-source federated learning platform that connects NHS Trusts and serves AI development across the UK population. I work end-to-end: foundation-model research (multimodal vision-language models for radiology), platform engineering (FastAPI, Docker Swarm, AWS, OpenTofu) and federated MLOps (NVIDIA FLARE, Flower).',
  metrics: [
    { value: '18', label: 'Q1 journal articles' },
    { value: '5,000+', label: 'Citations · h-index 19' }
  ],
  highlights: [
    {
      title: 'FLIP',
      description: 'Open-source federated learning platform connecting NHS Trusts.',
      href: 'https://github.com/londonaicentre/FLIP'
    },
    {
      title: 'Audify',
      description: 'LLM-powered Spotify playlist generator with semantic search.',
      href: 'https://github.com/garciadias/audify'
    },
    {
      title: 'Neurofind',
      description: 'Diagnostic tool aiding mental-disorder diagnosis from 3D brain MRI.',
      href: 'https://neurofind.ai/'
    },
    {
      title: 'Monzo overdraft GBM',
      description: '20–100% Gini gain on a £259M lending portfolio serving 7M+ customers.',
      href: 'https://monzo.com/static/docs/monzo-annual-report-2022.pdf'
    },
    {
      title: 'Neuroharmony',
      description: 'First out-of-sample harmonisation tool for neuroimaging.',
      href: 'https://neuroharmony.readthedocs.io/en/latest/'
    },
    {
      title: 'MONAI contributor',
      description: 'Active contributor to MONAI Core & MONAI Generative.',
      href: 'https://monai.io/'
    }
  ],
  stack: {
    'Languages & ML': ['Python', 'PyTorch', 'scikit-learn', 'TensorFlow', 'R', 'SQL'],
    'Foundation models / GenAI': ['Vision-language models', 'LLM fine-tuning & evaluation', 'RAG', 'OpenAI APIs', 'MONAI Generative'],
    'Federated learning': ['NVIDIA FLARE', 'Flower Framework', 'Transfer learning'],
    'Platforms & infra': ['AWS', 'Docker', 'Docker Swarm', 'Kubernetes', 'FastAPI', 'Django', 'OpenTofu/Terraform', 'GCP', 'BigQuery', 'dbt'],
    MLOps: ['CI/CD (GitHub Actions)', 'MLflow', 'Model monitoring', 'Experiment tracking'],
    'Healthcare data': ['DICOM', 'OMOP CDM', 'XNAT', 'Orthanc PACS', 'NHS information governance']
  },
  languages: ['English', 'Spanish', 'Portuguese']
}
