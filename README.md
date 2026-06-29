# BBB Peptide Score

A production-ready MVP Next.js portal for rapid **in silico prioritization of BBB-penetrating peptide candidates**. Users can analyze one peptide or a batch of named sequences, inspect transparent physicochemical metrics, review warnings and interpretation, and export results to CSV.

## Features

- Next.js App Router + TypeScript + Tailwind CSS
- Full-stack API route at `POST /api/analyze`
- Single and batch peptide input modes
- Validation for standard amino acid sequences
- Transparent metrics: length, molecular weight, net charge at pH 7.4, residue fractions, GRAVY, estimated pI, CPP-like score, solubility and toxicity-risk labels
- Heuristic BBB Penetration Score from 0–100 with Low / Medium / High class
- Confidence labels, warnings, interpretation bullets, validation suggestions
- Client-side CSV export
- No database, authentication, or external API calls
- Deployable on Vercel

## Install

```bash
npm install
```

## Run locally

```bash
npm run dev
```

Open <http://localhost:3000>.

## Build

```bash
npm run build
```

## Example test sequences

```text
TAT,YGRKKRRQRRR
Angiopep-2,TFFYGGSRGKRNNFKTEEY
PolyK,KKKKKKKKKK
AcidicControl,DEDEDEDEDE
HydrophobicControl,LLLLLLVVVVVV
```

## Scoring method summary

The score is a transparent heuristic for research triage. It starts at 50 and adjusts for:

- Net positive charge
- Arg/Lys/basic residue enrichment
- Length suitability, especially 6–30 amino acids
- Aromatic residue fraction
- Hydropathy balance using Kyte-Doolittle GRAVY
- Solubility-related indicators
- Penalties for short/long peptides, acidic enrichment, cysteine complexity, extreme cationic density, and extreme hydrophobicity

Approximate molecular weight uses average residue masses minus water loss per peptide bond. Net charge at pH 7.4 is a simple approximation using termini, K/R/H, and D/E contributions. Estimated pI is a rough proxy suitable only for MVP-level prioritization.

## Disclaimer

This tool provides an in silico prioritization score for research use only. It does not prove BBB penetration, brain exposure, pharmacokinetics, safety, or therapeutic efficacy. Experimental validation is required.
