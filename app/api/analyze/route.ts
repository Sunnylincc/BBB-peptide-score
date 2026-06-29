import { NextResponse } from 'next/server';
import { analyzePeptide, parseInput, type Mode } from '@/lib/peptide';

export async function POST(request: Request) {
  try {
    const body = (await request.json()) as { input?: string; mode?: Mode };
    const input = body.input ?? '';
    const mode = body.mode === 'batch' ? 'batch' : 'single';
    if (!input.trim()) return NextResponse.json({ error: 'Input is required.' }, { status: 400 });
    const entries = parseInput(input, mode);
    if (entries.length === 0) return NextResponse.json({ error: 'No peptide sequences were found.' }, { status: 400 });
    const results = entries.map((entry) => analyzePeptide(entry.name, entry.sequence)).sort((a, b) => b.bbbScore - a.bbbScore);
    return NextResponse.json({ results });
  } catch {
    return NextResponse.json({ error: 'Unable to analyze peptide input.' }, { status: 500 });
  }
}
