import type { Metadata } from 'next';
import './globals.css';

export const metadata: Metadata = { title: 'BBB Peptide Score', description: 'Rapid in silico prioritization of BBB-penetrating peptide candidates' };
export default function RootLayout({ children }: { children: React.ReactNode }) { return <html lang="en"><body>{children}</body></html>; }
