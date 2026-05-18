import html2pdf from 'html2pdf.js'

export function generateCVPDF(experience, education, profile) {
  const htmlContent = `
    <div style="font-family: Arial, sans-serif; color: #000; line-height: 1.3;">
      <!-- Header -->
      <div style="margin-bottom: 10px; border-bottom: 2px solid #000; padding-bottom: 8px;">
        <h1 style="margin: 0 0 2px 0; font-size: 20px; font-weight: bold;">${profile.shortName}</h1>
        <p style="margin: 2px 0; font-size: 10px; color: #333;">
          ${profile.location} • ${profile.email}
        </p>
      </div>

      <!-- Professional Summary -->
      <div style="margin-bottom: 12px;">
        <h2 style="margin: 0 0 4px 0; font-size: 11px; font-weight: bold; text-transform: uppercase; letter-spacing: 0.5px;">Professional Summary</h2>
        <p style="margin: 0; font-size: 10px; line-height: 1.3; color: #333;">
          ${profile.intro}
        </p>
      </div>

      <!-- Experience -->
      <div style="margin-bottom: 12px;">
        <h2 style="margin: 0 0 8px 0; font-size: 11px; font-weight: bold; text-transform: uppercase; letter-spacing: 0.5px;">Professional Experience</h2>
        ${experience.map((role, idx) => `
          <div style="margin-bottom: 10px; ${idx > 0 ? 'border-top: 1px solid #ddd; padding-top: 8px;' : ''}">
            <div style="display: flex; justify-content: space-between; align-items: baseline; margin-bottom: 1px;">
              <h3 style="margin: 0; font-size: 10px; font-weight: bold;">${role.role}</h3>
              <span style="font-size: 9px; color: #666; margin-left: 12px; white-space: nowrap;">${role.period}</span>
            </div>
            <p style="margin: 1px 0 3px 0; font-size: 9px; color: #555;">
              ${role.org}${role.location ? ' • ' + role.location : ''}
            </p>
            ${role.summary ? `<p style="margin: 2px 0 3px 0; font-size: 9px; color: #333;">${role.summary}</p>` : ''}
            ${role.sections && role.sections.length > 0 ? `
              <ul style="margin: 3px 0 0 0; padding-left: 14px;">
                ${role.sections.flatMap(s => 
                  s.bullets ? s.bullets.map(bullet => `
                    <li style="margin: 1px 0; font-size: 9px; color: #333; line-height: 1.3;">${bullet}</li>
                  `).join('') : ''
                ).join('')}
              </ul>
            ` : ''}
          </div>
        `).join('')}
      </div>

      <!-- Education -->
      <div style="margin-bottom: 12px;">
        <h2 style="margin: 0 0 6px 0; font-size: 11px; font-weight: bold; text-transform: uppercase; letter-spacing: 0.5px;">Education</h2>
        ${education.map((edu, idx) => `
          <div style="margin-bottom: 4px; ${idx > 0 ? 'border-top: 1px solid #ddd; padding-top: 4px;' : ''}">
            <div style="display: flex; justify-content: space-between; align-items: baseline;">
              <h3 style="margin: 0; font-size: 10px; font-weight: bold;">${edu.degree}</h3>
              <span style="font-size: 9px; color: #666; margin-left: 12px;">${edu.year}</span>
            </div>
            <p style="margin: 1px 0 0 0; font-size: 9px; color: #555;">${edu.institution}</p>
          </div>
        `).join('')}
      </div>

      <!-- Languages -->
      ${profile.languages && profile.languages.length > 0 ? `
        <div>
          <h2 style="margin: 0 0 4px 0; font-size: 11px; font-weight: bold; text-transform: uppercase; letter-spacing: 0.5px;">Languages</h2>
          <p style="margin: 0; font-size: 9px; color: #333;">${profile.languages.join(' • ')}</p>
        </div>
      ` : ''}
    </div>
  `

  const opt = {
    margin: [8, 8, 8, 8],
    filename: 'CV-Rafael-Garcia-Dias.pdf',
    image: { type: 'jpeg', quality: 0.98 },
    html2canvas: { scale: 2, useCORS: true },
    jsPDF: { unit: 'mm', format: 'a4', orientation: 'portrait' },
    pagebreak: { mode: ['avoid-all', 'css', 'legacy'] }
  }

  html2pdf().set(opt).from(htmlContent).save()
}
