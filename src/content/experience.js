export const experience = [
  {
    role: 'Staff AI Engineer — Foundation Models for Healthcare',
    org: 'AI Centre for Value-Based Healthcare / King’s College London',
    location: 'London, UK',
    period: 'Nov 2024 – Present',
    summary:
      'Principal technical lead and architect for production AI infrastructure and foundation models across an NHS Trust federation, in collaboration with universities (KCL, GSTT) and tier-1 industry partners (NVIDIA, deepc, Flower Labs, OneLondon).',
    sections: [
      {
        heading: 'FLIP — Federated Learning & Interoperability Platform',
        subheading: 'Open source · Apache 2.0',
        href: 'https://github.com/londonaicentre/FLIP',
        bullets: [
          'Platform architecture & strategy: architected and served as the principal technical lead for FLIP, a multi-Trust federated learning platform (publicly launched Feb 2026), defining the technical strategy for privacy-preserving AI training across NHS firewalls without moving patient data.',
          'Cross-functional leadership: acted as the primary technical liaison aligning clinical researchers, academic institutions (KCL, GSTT) and tier-1 industry partners (NVIDIA, Flower Labs) to establish standardised, secure deployment workflows.',
          'System design: designed a highly scalable microservice architecture (FastAPI, Docker Swarm, AWS) featuring a central hub API, per-Trust APIs, and secure outbound-only Trust-to-Hub communication via CloudFront/ALB; contributed the OpenTofu/Terraform AWS deployment for production environments.',
          'MLOps & framework integration: directed the integration of dual industry-standard federated learning backends (NVIDIA FLARE and Flower Framework), enabling parallel execution and standardising the data plane for diverse research teams.',
          'Integrated DICOM imaging infrastructure (Orthanc, XNAT) and OMOP CDM structured-data layer, giving researchers unified multimodal access across imaging and EHR within each Trust’s secure enclave.'
        ]
      },
      {
        heading: 'Multimodal foundation models for radiology',
        bullets: [
          'Foundational model innovation: spearheaded the development of cutting-edge vision-language foundation models fusing multi-modal medical imaging (X-ray, CT, MRI) with radiology reports — supporting automated report generation, image synthesis from reports, scan/report discrepancy detection, and LLM-based report classification — setting the methodology for responsible, clinically safe LLM scaling.',
          'Validated robustness of latent representations (including LLM embeddings) so smaller task-specific models built on top remain reliable — a prerequisite for clinical translation.',
          'Active open-source contributor to MONAI Core and MONAI Generative, advancing foundational algorithms used by the wider medical-imaging research community.',
          'Drives methodology for the responsible scaling of LLMs in clinical settings — bridging cutting-edge AI capability with NHS governance, safety and real-world clinical workflows.'
        ]
      }
    ]
  },
  {
    role: 'Founding ML Engineer',
    org: 'Floe Oral Care',
    location: 'London, UK',
    period: 'May 2023 – Nov 2024',
    summary: 'Founding ML engineer; defined and executed the end-to-end technical architecture, taking the product from research concept to deployed clinical diagnostic tool.',
    sections: [
      {
        bullets: [
          'Technical vision & execution: defined the end-to-end technical architecture, taking the product from an initial research concept to a fully deployed clinical diagnostic tool.',
          'Infrastructure & standards: architected the complete cloud stack on AWS (CDK), establishing rigorous CI/CD pipelines, SOLID principles, and a strong organisational testing culture.',
          'Engineering leadership: directed the front-end engineering team and built the integration layers (Django backend) to ensure seamless communication between user-facing applications and ML assets.',
          'Owned the science: designed the clinical trials that validated the diagnostic concept and produced the training data for the production models.',
          'Shipped multiple ML models: logistic regression mapping protein concentrations to health risk, tree-based models over questionnaire data, and an OpenAI-backed RAG recommendation system; iterated production models post-launch based on real user data.'
        ]
      }
    ]
  },
  {
    role: 'Decision Scientist (Credit Risk ML)',
    org: 'Monzo Bank',
    location: 'London, UK',
    period: 'Aug 2022 – Feb 2023 · Fixed-term, 6 months',
    summary:
      'Improved creditworthiness modelling on a £259M lending portfolio serving 7M+ customers.',
    sections: [
      {
        bullets: [
          'Organisational standardisation: spearheaded the refactoring of the internal modelling library to standardise ML monitoring practices across all Monzo lending products, surfacing critical performance gaps across previously unmonitored customer subgroups.',
          'High-impact delivery: shipped a GBM-based overdraft model (GCP/Docker) that achieved 20–100% Gini gains over previous baselines on a £259M lending portfolio serving 7M+ customers.',
          'Stack: Python, scikit-learn, BigQuery, dbt, Looker, GCP.'
        ]
      }
    ]
  },
  {
    role: 'Research Associate — Founding Engineer, Neurofind',
    org: 'Neurofind / King’s College London',
    location: 'London, UK',
    period: 'Aug 2018 – Aug 2022 · 4 years, fixed-term renewed 4×',
    summary: '',
    sections: [
      {
        bullets: [
          'Built the backend, cloud infrastructure and ML models for Neurofind.ai, a tool aiding diagnosis of mental disorders from 3D brain MRI.',
          'Authored Neuroharmony, the first published out-of-sample harmonisation tool for neuroimaging — a key step toward clinical translation of imaging ML by mitigating cross-scanner bias.',
          'Won a £109,000 MRC grant as named investigator.',
          'Authored one chapter and co-authored five more in Machine Learning: Methods and Applications to Brain Disorders (Elsevier).',
          'Stack: Python, PyTorch, TensorFlow, scikit-learn, AWS, Docker, MLflow, Plotly, R.'
        ]
      }
    ]
  },
  {
    role: 'PhD Fellow — Machine Learning Applied to Astrophysics',
    org: 'Instituto de Astrofísica de Canarias (IAC) / Universidad de La Laguna',
    location: 'Tenerife, Spain',
    period: 'Aug 2015 – Aug 2018',
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
