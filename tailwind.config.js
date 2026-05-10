/** @type {import('tailwindcss').Config} */
export default {
  content: ['./index.html', './src/**/*.{vue,js,ts,jsx,tsx}'],
  darkMode: 'class',
  theme: {
    extend: {
      fontFamily: {
        sans: ['Inter', 'ui-sans-serif', 'system-ui', 'sans-serif'],
        mono: ['JetBrainsMono', 'ui-monospace', 'SFMono-Regular', 'monospace'],
        display: ['"Bai Jamjuree"', 'Inter', 'sans-serif']
      },
      colors: {
        // Primary accent — Dracula purple (#bd93f9) shaded for full Tailwind range
        primary: {
          50: '#faf6ff',
          100: '#f3eaff',
          200: '#e8d6ff',
          300: '#d6b6fc',
          400: '#bd93f9',
          500: '#a370f0',
          600: '#8a4ee0',
          700: '#743cc4',
          800: '#5e329e',
          900: '#4a2a7a',
          950: '#2d1850'
        },
        // Dracula palette — for surfaces, accents, and code
        dracula: {
          bg: '#282a36',
          surface: '#21222c',
          line: '#44475a',
          fg: '#f8f8f2',
          comment: '#6272a4',
          pink: '#ff79c6',
          purple: '#bd93f9',
          cyan: '#8be9fd',
          green: '#50fa7b',
          yellow: '#f1fa8c',
          orange: '#ffb86c',
          red: '#ff5555'
        }
      },
      typography: {
        DEFAULT: {
          css: {
            maxWidth: 'none'
          }
        }
      }
    }
  },
  plugins: []
}
