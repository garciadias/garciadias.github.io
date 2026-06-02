export const experience = [
  {
    role: 'Senior AI Engineer - Foundation Models for Healthcare',
    org: 'AI Centre for Value-Based Healthcare / King’s College London',
    location: 'London, UK',
    period: 'Nov 2024 - Present',
    summary:
      'Building production AI infrastructure, distributed training systems, and foundation models across three large NHS Trusts encompassing millions of patients, in collaboration with KCL, GSTT, NVIDIA, deepc, Flower Labs, and OneLondon.',
    sections: [
      {
        heading: 'FLIP - Federated Learning & Interoperability Platform',
        subheading: 'Open source · Apache 2.0',
        href: 'https://github.com/londonaicentre/FLIP',
        bullets: [
          'Core contributor to the multi-Trust platform launched publicly in Feb 2026, enabling privacy-preserving distributed training of AI across three NHS Trusts (GSTT, KCH, BDMS) without moving patient data. Microservice architecture: central hub API, per-Trust APIs, imaging and data-access services, FastAPI throughout, deployed via Docker Swarm and AWS.',
          'Appointed Federated Learning Chair of the MONAI Working Group, co-chairing alongside NVIDIA to set direction for federated learning standards and best practices across the global medical-imaging community.',
          'Integrated NVIDIA FLARE and Flower Framework as parallel federated-learning backends, allowing research teams to choose between two industry-standard frameworks against the same data plane.',
          'Designed the trust-authentication and internal-service-key model for secure outbound-only Trust-to-Hub communication via CloudFront/ALB; contributed to the OpenTofu/Terraform AWS deployment for production environments.',
          'Integrated DICOM imaging infrastructure (Orthanc, XNAT) and OMOP CDM structured-data layer, giving researchers unified multimodal access across imaging and EHR within each Trust’s secure enclave.'
        ]
      },
      {
        heading: 'Multimodal foundation models for radiology',
        bullets: [
          'Leading development of vision-language foundation models that fuse medical imaging (X-ray, CT, MRI) with radiology reports, supporting downstream applications: automated report generation, image synthesis from reports, scan/report discrepancy detection, and LLM-based report classification.',
          'Validated robustness of latent representations (including LLM embeddings) so smaller task-specific models built on top remain reliable - a prerequisite for clinical translation.',
          'Active open-source contributor to MONAI Core and MONAI Generative, advancing foundational algorithms used by the wider medical-imaging research community.',
          'Drives methodology for the responsible scaling of LLMs in clinical settings — bridging cutting-edge AI capability with NHS governance, safety and real-world clinical workflows.',
          'Applied transfer learning by fine-tuning foundation models on the myskin.org psoriasis detection task, demonstrating how models pre-trained on broad biomedical data adapt to a specific downstream clinical diagnostic use case.'
        ]
      }
    ]
  },
  {
    role: 'Machine Learning Engineer',
    org: 'Floe Oral Care',
    location: 'London, UK',
    period: 'May 2023 - Nov 2024',
    summary: 'Founding ML engineer; took the product from research concept to deployed clinical tool.',
    sections: [
      {
        bullets: [
          'Owned the science: designed the clinical trials that validated the diagnostic concept and produced the training data for the production models.',
          'Shipped multiple ML models: logistic regression mapping protein concentrations to health risk, tree-based models over questionnaire data, and an OpenAI-backed RAG recommendation system.',
          'Built the full cloud stack on AWS with AWS CDK, including databases, ML asset storage and web infrastructure; CI/CD throughout.',
          'Created the Django backend, integrated the contractor-built front end and directed the front-end engineering team.',
          'Iterated production models post-launch based on real user data, applying SOLID principles and a strong testing culture.'
        ]
      }
    ]
  },
  {
    role: 'Decision Scientist (Credit Risk ML)',
    org: 'Monzo Bank',
    location: 'London, UK',
    period: 'Aug 2022 - Feb 2023 · Fixed-term, 6 months',
    summary:
      'Improved creditworthiness modelling on a £259M lending portfolio serving 7M+ customers.',
    sections: [
      {
        bullets: [
          'Shipped a GBM-based overdraft-application model delivering 20-100% Gini gains over the previous baseline, deployed in Docker on GCP.',
          'Refactored the internal modelling library to standardise monitoring across Monzo’s lending products, surfacing performance gaps across customer subgroups that had previously been invisible.',
          'Mentored a junior team member, guiding their development in production ML engineering and credit-risk modelling.',
          'Stack: Python, scikit-learn, BigQuery, dbt, Looker, GCP.'
        ]
      }
    ]
  },
  {
    role: 'Research Associate - Founding Engineer, Neurofind',
    org: 'Neurofind / King’s College London',
    location: 'London, UK',
    period: 'Aug 2018 - Aug 2022 · 4 years, fixed-term renewed 4×',
    summary: '',
    sections: [
      {
        bullets: [
          'Built the backend, cloud infrastructure and ML models for Neurofind.ai, a tool aiding diagnosis of mental disorders from 3D brain MRI.',
          'Authored Neuroharmony, the first published out-of-sample harmonisation tool for neuroimaging - a key step toward clinical translation of imaging ML by mitigating cross-scanner bias.',
          'Won a £109,000 MRC grant as named investigator.',
          'Authored one chapter and co-authored five more in Machine Learning: Methods and Applications to Brain Disorders (Elsevier).',
          'Supervised a Master’s and a PhD student, guiding their research in neuroimaging ML and contributing to the group’s publication output.',
          'Stack: Python, PyTorch, TensorFlow, scikit-learn, AWS, Docker, MLflow, Plotly, R.'
        ]
      }
    ]
  },
  {
    role: 'PhD Fellow - Machine Learning Applied to Astrophysics',
    org: 'Instituto de Astrofísica de Canarias (IAC) / Universidad de La Laguna',
    location: 'Tenerife, Spain',
    period: 'Aug 2015 - Aug 2018',
    summary:
      'Applied unsupervised and supervised ML to large-scale stellar spectra surveys; developed clustering methodology that became the basis for a later book chapter and journal paper. Foundation in rigorous, reproducible scientific computing in Python (NumPy, pandas, scikit-learn, Matplotlib) on Linux.',
    sections: []
  }
]

export const education = [
  {
    degree: 'PhD, Machine Learning Applied to Astrophysics',
    institution: 'Instituto de Astrofísica de Canarias / Universidad de La Laguna',
    year: '2018'
  }
]
