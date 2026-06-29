export type Mode = 'single' | 'batch';
export type PredictionClass = 'Low' | 'Medium' | 'High';
export type Confidence = 'Low' | 'Medium' | 'High';

export interface PeptideMetrics {
  length: number;
  molecularWeight: number;
  netChargePh74: number;
  basicFraction: number;
  acidicFraction: number;
  hydrophobicFraction: number;
  aromaticFraction: number;
  polarFraction: number;
  gravy: number;
  estimatedPI: number;
  cppLikeScore: number;
  solubilityLabel: string;
  toxicityRiskLabel: string;
}

export interface PeptideResult {
  name: string;
  sequence: string;
  valid: boolean;
  bbbScore: number;
  predictionClass: PredictionClass;
  confidence: Confidence;
  metrics: PeptideMetrics;
  warnings: string[];
  interpretation: string[];
  validationSuggestions: string[];
}

const ALLOWED_SEQUENCE = 'ACDEFGHIKLMNPQRSTVWY';
const ALLOWED_AMINO_ACIDS = new Set(ALLOWED_SEQUENCE.split(''));
const AVERAGE_RESIDUE_MASSES: Record<string, number> = {
  A: 89.09, R: 174.20, N: 132.12, D: 133.10, C: 121.15,
  E: 147.13, Q: 146.15, G: 75.07, H: 155.16, I: 131.17,
  L: 131.17, K: 146.19, M: 149.21, F: 165.19, P: 115.13,
  S: 105.09, T: 119.12, W: 204.23, Y: 181.19, V: 117.15,
};
const KYTE_DOOLITTLE: Record<string, number> = {
  A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5,
  Q: -3.5, E: -3.5, G: -0.4, H: -3.2, I: 4.5,
  L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6,
  S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2,
};

const HYDROPHOBIC = new Set('AVILMFWYP'.split(''));
const AROMATIC = new Set('FWY'.split(''));
const BASIC = new Set('KRH'.split(''));
const ACIDIC = new Set('DE'.split(''));
const POLAR = new Set('STNQCY'.split(''));
const CYSTEINE = new Set(['C']);

export const validationSuggestions = [
  'Brain endothelial cell uptake assay',
  'In vitro BBB transwell assay',
  'Serum stability assay',
  'Cytotoxicity or hemolysis assay',
  'In vivo biodistribution if early assays pass',
];

export function parseInput(input: string, mode: Mode): { name: string; sequence: string }[] {
  const lines = input.split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
  const selectedLines = mode === 'single' ? lines.slice(0, 1) : lines;

  return selectedLines.map((line, index) => {
    const [first, ...rest] = line.split(',');
    if (rest.length === 0) {
      return { name: `peptide_${index + 1}`, sequence: first ?? '' };
    }

    return {
      name: (first || `peptide_${index + 1}`).trim(),
      sequence: rest.join(','),
    };
  });
}

export function cleanSequence(sequence: string): string {
  return sequence.replace(/\s+/g, '').toUpperCase();
}

export function validateSequence(sequence: string): {
  sequence: string;
  valid: boolean;
  hasNonStandard: boolean;
} {
  const cleanedSequence = cleanSequence(sequence);
  const hasNonStandard = [...cleanedSequence].some((aa) => !ALLOWED_AMINO_ACIDS.has(aa));

  return {
    sequence: cleanedSequence,
    valid: cleanedSequence.length > 0 && !hasNonStandard,
    hasNonStandard,
  };
}

function residueFraction(sequence: string, residues: Set<string>): number {
  return sequence.length ? [...sequence].filter((aa) => residues.has(aa)).length / sequence.length : 0;
}

function round(value: number, digits: number): number {
  return Number(value.toFixed(digits));
}

export function calculateFeatures(rawSequence: string): PeptideMetrics {
  const sequence = cleanSequence(rawSequence).replace(/[^ACDEFGHIKLMNPQRSTVWY]/g, '');
  const length = sequence.length;
  const massSum = [...sequence].reduce((sum, aa) => sum + AVERAGE_RESIDUE_MASSES[aa], 0);

  // Approximate peptide MW using average residue masses minus water for each peptide bond.
  const molecularWeight = length ? massSum - 18.015 * (length - 1) : 0;

  // Approximate net charge at pH 7.4: termini plus K/R, partial H, and D/E contributions.
  const netChargePh74 = length
    ? [...sequence].reduce((sum, aa) => {
        if (aa === 'K' || aa === 'R') return sum + 1;
        if (aa === 'H') return sum + 0.1;
        if (aa === 'D' || aa === 'E') return sum - 1;
        return sum;
      }, 0)
    : 0;

  const basicFraction = residueFraction(sequence, BASIC);
  const acidicFraction = residueFraction(sequence, ACIDIC);
  const hydrophobicFraction = residueFraction(sequence, HYDROPHOBIC);
  const aromaticFraction = residueFraction(sequence, AROMATIC);
  const polarFraction = residueFraction(sequence, POLAR);
  const gravy = length ? [...sequence].reduce((sum, aa) => sum + KYTE_DOOLITTLE[aa], 0) / length : 0;

  const estimatedPI = estimatePI(netChargePh74, basicFraction);
  const cppLikeScore = calculateCppLikeScore({ length, netChargePh74, basicFraction, gravy, hydrophobicFraction });

  return {
    length,
    molecularWeight: round(molecularWeight, 1),
    netChargePh74: round(netChargePh74, 1),
    basicFraction: round(basicFraction, 2),
    acidicFraction: round(acidicFraction, 2),
    hydrophobicFraction: round(hydrophobicFraction, 2),
    aromaticFraction: round(aromaticFraction, 2),
    polarFraction: round(polarFraction, 2),
    gravy: round(gravy, 2),
    estimatedPI: round(estimatedPI, 1),
    cppLikeScore: round(cppLikeScore, 2),
    solubilityLabel: getSolubilityLabel(hydrophobicFraction, gravy, netChargePh74),
    toxicityRiskLabel: getToxicityRiskLabel(hydrophobicFraction, netChargePh74, basicFraction, length),
  };
}

function estimatePI(netChargePh74: number, basicFraction: number): number {
  if (netChargePh74 > 4 && basicFraction > 0.25) return Math.min(12.5, 10 + basicFraction * 3);
  if (netChargePh74 < -2) return Math.max(3.5, 5.5 + netChargePh74 * 0.2);
  return 7 + Math.max(-1, Math.min(1, netChargePh74 * 0.15));
}

function calculateCppLikeScore(metrics: {
  length: number;
  netChargePh74: number;
  basicFraction: number;
  gravy: number;
  hydrophobicFraction: number;
}): number {
  const chargeScore = Math.min(1, Math.max(0, metrics.netChargePh74 / 8));
  const basicScore = Math.min(1, metrics.basicFraction / 0.5);
  const lengthScore = metrics.length >= 6 && metrics.length <= 30
    ? 1
    : metrics.length < 6
      ? metrics.length / 6
      : Math.max(0, 1 - (metrics.length - 30) / 30);
  const hydropathyScore = metrics.gravy <= 2 && metrics.hydrophobicFraction < 0.6 ? 1 : 0.35;

  return chargeScore * 0.35 + basicScore * 0.3 + lengthScore * 0.2 + hydropathyScore * 0.15;
}

function getSolubilityLabel(hydrophobicFraction: number, gravy: number, netChargePh74: number): string {
  if (hydrophobicFraction > 0.6 || gravy > 2) return 'Aggregation risk';
  if (netChargePh74 > 0 && hydrophobicFraction < 0.45) return 'Likely soluble';
  return 'Moderate solubility';
}

function getToxicityRiskLabel(
  hydrophobicFraction: number,
  netChargePh74: number,
  basicFraction: number,
  length: number,
): string {
  if (netChargePh74 / Math.max(1, length) > 0.75 || basicFraction > 0.75) return 'High cationic risk';
  if (netChargePh74 > 4) return 'Moderate cationic risk';
  if (hydrophobicFraction > 0.65) return 'Hydrophobic membrane risk';
  return 'Lower heuristic risk';
}

export function generateWarnings(sequence: string, metrics: PeptideMetrics, hasNonStandard = false): string[] {
  const warnings: string[] = [];
  const cysteineFraction = residueFraction(sequence, CYSTEINE);

  if (hasNonStandard) warnings.push('Sequence contains non-standard amino acid characters.');
  if (metrics.length < 5) warnings.push('Very short peptides may be unstable or nonspecific.');
  if (metrics.length > 50) warnings.push('Long peptides may have reduced permeability and higher synthesis complexity.');
  if (metrics.gravy > 2 || metrics.hydrophobicFraction > 0.65) warnings.push('Highly hydrophobic peptides may have poor solubility or aggregation risk.');
  if (metrics.acidicFraction > 0.3) warnings.push('Highly acidic peptides are deprioritized by this BBB/CPP heuristic.');
  if (metrics.netChargePh74 / Math.max(1, metrics.length) > 0.75 || metrics.basicFraction > 0.75 || metrics.netChargePh74 > 8) warnings.push('Extremely cationic peptides may show nonspecific binding or cytotoxicity.');
  if (cysteineFraction > 0.15) warnings.push('High cysteine content may cause oxidation or disulfide-related complexity.');

  return [...new Set(warnings)];
}

export function calculateBBBScore(metrics: PeptideMetrics): number {
  let score = 50;
  const charge = metrics.netChargePh74;

  if (charge >= 1 && charge <= 3) score += 6;
  else if (charge >= 4 && charge <= 8) score += 14;
  else if (charge > 8) score += 16;

  if (metrics.basicFraction >= 0.15 && metrics.basicFraction <= 0.45) score += 8;
  else if (metrics.basicFraction > 0.45 && metrics.basicFraction <= 0.75) score += 15;
  else if (metrics.basicFraction > 0.75) score += 10;

  if (metrics.length >= 6 && metrics.length <= 20) score += 15;
  else if (metrics.length >= 21 && metrics.length <= 35) score += 10;
  else if (metrics.length >= 36 && metrics.length <= 50) score += 5;

  if (metrics.aromaticFraction >= 0.05 && metrics.aromaticFraction <= 0.25) score += 8;
  else if (metrics.aromaticFraction > 0.25) score += 5;

  if (metrics.gravy >= -2.5 && metrics.gravy <= 1) score += 10;
  else if (metrics.gravy > 1 && metrics.gravy <= 2) score += 5;
  else if (metrics.gravy > 2) score -= 10;
  else if (metrics.gravy < -3) score -= 5;

  if (metrics.hydrophobicFraction < 0.45 && metrics.netChargePh74 > 0) score += 10;
  else if (metrics.hydrophobicFraction > 0.6) score -= 10;

  if (metrics.length < 5) score -= 15;
  if (metrics.length > 50) score -= 15;
  if (metrics.acidicFraction > 0.3) score -= 15;
  if (metrics.netChargePh74 / Math.max(1, metrics.length) > 0.75) score -= 5;
  if (metrics.hydrophobicFraction > 0.65) score -= 15;

  return Math.max(0, Math.min(100, Math.round(score)));
}

export function predictionClass(score: number): PredictionClass {
  if (score < 40) return 'Low';
  if (score < 70) return 'Medium';
  return 'High';
}

export function confidence(valid: boolean, metrics: PeptideMetrics, warnings: string[], score: number): Confidence {
  const nearBoundary = [39, 40, 69, 70].some((boundary) => Math.abs(score - boundary) <= 4);
  if (valid && metrics.length >= 6 && metrics.length <= 30 && warnings.length === 0 && !nearBoundary) return 'High';
  if (!valid || warnings.length > 0 || metrics.length < 5 || metrics.length > 50 || metrics.gravy > 2 || metrics.basicFraction > 0.75) return 'Low';
  return 'Medium';
}

export function generateInterpretation(metrics: PeptideMetrics): string[] {
  const interpretation = [
    metrics.netChargePh74 > 4
      ? 'Strongly cationic sequence, consistent with CPP-like uptake behavior but requiring toxicity checks.'
      : metrics.netChargePh74 > 0
        ? 'Moderately positive charge may support interaction with anionic cell surfaces.'
        : 'Neutral or negative charge lowers this heuristic BBB/CPP prioritization score.',
    metrics.length >= 6 && metrics.length <= 20
      ? 'Short peptide length is compatible with BBB shuttle-like candidates.'
      : 'Peptide length may increase permeability or synthesis uncertainty.',
    metrics.gravy < -2.5
      ? 'Very hydrophilic profile supports solubility but may reduce passive membrane diffusion.'
      : metrics.gravy > 2
        ? 'Hydrophobic profile may increase aggregation or nonspecific membrane interactions.'
        : 'Balanced hydropathy is favorable for early prioritization.',
  ];

  if (metrics.aromaticFraction >= 0.05) {
    interpretation.push('Aromatic residues may contribute membrane interaction and receptor-binding opportunities.');
  }

  interpretation.push(`CPP-like score is ${metrics.cppLikeScore.toFixed(2)}, a transparent heuristic rather than experimental evidence.`);
  return interpretation.slice(0, 5);
}

export function analyzePeptide(name: string, rawSequence: string): PeptideResult {
  const validation = validateSequence(rawSequence);
  const metrics = calculateFeatures(validation.sequence);
  const warnings = generateWarnings(validation.sequence, metrics, validation.hasNonStandard);
  const bbbScore = validation.valid ? calculateBBBScore(metrics) : 0;

  return {
    name,
    sequence: validation.sequence,
    valid: validation.valid,
    bbbScore,
    predictionClass: predictionClass(bbbScore),
    confidence: confidence(validation.valid, metrics, warnings, bbbScore),
    metrics,
    warnings,
    interpretation: generateInterpretation(metrics),
    validationSuggestions,
  };
}

export function generateCSV(results: PeptideResult[]): string {
  const headers = [
    'name', 'sequence', 'valid', 'bbbScore', 'predictionClass', 'confidence', 'length',
    'molecularWeight', 'netChargePh74', 'basicFraction', 'hydrophobicFraction',
    'aromaticFraction', 'gravy', 'estimatedPI', 'cppLikeScore', 'warnings',
  ];
  const rows = results.map((result) => [
    result.name,
    result.sequence,
    result.valid,
    result.bbbScore,
    result.predictionClass,
    result.confidence,
    result.metrics.length,
    result.metrics.molecularWeight,
    result.metrics.netChargePh74,
    result.metrics.basicFraction,
    result.metrics.hydrophobicFraction,
    result.metrics.aromaticFraction,
    result.metrics.gravy,
    result.metrics.estimatedPI,
    result.metrics.cppLikeScore,
    result.warnings.join('; '),
  ]);

  return [headers, ...rows]
    .map((row) => row.map((cell) => `"${String(cell).replace(/"/g, '""')}"`).join(','))
    .join('\n');
}
