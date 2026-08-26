const fs = require("fs");
const path = require("path");
const D = require("docx");

const blocks = JSON.parse(fs.readFileSync("blocks.json", "utf8"));
const FIGDIR = "D:/R_ex/MR/figures";
const FIGS = [
  ["Fig1_locus_attribution", "Locus attribution under both outcomes, by bounded locus."],
  ["Fig2_generality", "Generality: five nested releases, the transfer test, the second disease, the non-cancer outcome, the second exposure resource, and the crossed grid with its mismatched-locus control."],
  ["Fig3_power_stability", "Power and list stability by locus class, and the same stratification against full-power |z|."],
  ["Fig4_coloc_heidi", "Colocalisation versus SMR/HEIDI, with in-sample fine-mapping."],
  ["Fig5_instrument_ladder", "Instrument availability across the pathway at three levels."],
  ["Fig6_tpi1_window", "The activation window in which TPI1 is instrumentable: effect size against precision across the eight profiles."],
  ["Fig7_axis_chromatin", "The CD4\u207A metabolic axis and its chromatin signature."],
  ["Fig8_compartment", "Compartment attribution."],
  ["Fig9_patients", "Patients across three cohorts."],
];

// --- inline **bold** / *italic* / `code` -> TextRun[] ------------------------
function runs(text, base = {}) {
  const out = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) out.push(new D.TextRun({ ...base, text: text.slice(last, m.index) }));
    const tok = m[0];
    if (tok.startsWith("**")) out.push(new D.TextRun({ ...base, text: tok.slice(2, -2), bold: true }));
    else if (tok.startsWith("`")) out.push(new D.TextRun({ ...base, text: tok.slice(1, -1), font: "Consolas" }));
    else out.push(new D.TextRun({ ...base, text: tok.slice(1, -1), italics: true }));
    last = m.index + tok.length;
  }
  if (last < text.length) out.push(new D.TextRun({ ...base, text: text.slice(last) }));
  return out.length ? out : [new D.TextRun({ ...base, text: "" })];
}

const DOUBLE = { line: 480, after: 0 };          // 240 = single
const H = [null, D.HeadingLevel.HEADING_1, D.HeadingLevel.HEADING_1,
           D.HeadingLevel.HEADING_2, D.HeadingLevel.HEADING_3];

// The manuscript now carries a full legend per figure in its own Figures
// section. Emit those with the images at the end rather than twice: once as
// running text and once as a hardcoded one-liner under the plate.
const figStart = blocks.findIndex(x => x.t === "h" && /^Figures$/.test(x.text || ""));
const figEnd = blocks.findIndex((x, i) => i > figStart && x.t === "h" && /^Methods$/.test(x.text || ""));
const LEGEND = {};
const FIGNOTE = [];
if (figStart >= 0 && figEnd > figStart) {
  for (let i = figStart + 1; i < figEnd; i++) {
    const x = blocks[i];
    if (x.t !== "p") continue;
    const m = /^\*\*Fig\.\s*(\d+)\s*\|/.exec(x.text);
    if (m) LEGEND[Number(m[1])] = x.text;
    else FIGNOTE.push(x.text);
  }
}

const children = [];
let title = null;

for (let bi = 0; bi < blocks.length; bi++) {
  const b = blocks[bi];
  if (figStart >= 0 && bi >= figStart && bi < figEnd) continue;
  if (b.t === "h" && b.level === 1 && title === null) { title = b.text; continue; }

  if (b.t === "h") {
    children.push(new D.Paragraph({
      children: runs(b.text, { bold: true }),
      heading: H[b.level],
      spacing: { before: 320, after: 160 },
      keepNext: true,
    }));
  } else if (b.t === "p") {
    children.push(new D.Paragraph({
      children: runs(b.text),
      spacing: DOUBLE,
      alignment: D.AlignmentType.LEFT,
    }));
  } else if (b.t === "ul") {
    for (const it of b.items) {
      children.push(new D.Paragraph({
        children: runs(it),
        bullet: { level: 0 },
        spacing: { line: 360, after: 0 },
      }));
    }
  } else if (b.t === "table") {
    const cols = Math.max(...b.rows.map(r => r.length));
    const total = 9360;                               // 6.5" in DXA
    const w = Math.floor(total / cols);
    const colWidths = Array(cols).fill(w);
    colWidths[cols - 1] = total - w * (cols - 1);
    children.push(new D.Table({
      columnWidths: colWidths,
      width: { size: total, type: D.WidthType.DXA },
      rows: b.rows.map((r, ri) => new D.TableRow({
        children: Array.from({ length: cols }, (_, ci) => new D.TableCell({
          width: { size: colWidths[ci], type: D.WidthType.DXA },
          shading: ri === 0
            ? { type: D.ShadingType.CLEAR, fill: "EFEFEF" }
            : undefined,
          margins: { top: 60, bottom: 60, left: 100, right: 100 },
          children: [new D.Paragraph({
            children: runs(r[ci] || "", { size: 20 }),
            spacing: { line: 240, after: 0 },
          })],
        })),
        tableHeader: ri === 0,
      })),
    }));
    children.push(new D.Paragraph({ text: "", spacing: { after: 160 } }));
  }
}

// --- figures, each on its own page at the end -------------------------------
children.push(new D.Paragraph({ children: [new D.PageBreak()] }));
children.push(new D.Paragraph({
  children: [new D.TextRun({ text: "Figures", bold: true })],
  heading: D.HeadingLevel.HEADING_1,
  spacing: { after: 240 },
}));

FIGS.forEach(([stem, caption], idx) => {
  const p = path.join(FIGDIR, stem + ".png");
  if (!fs.existsSync(p)) { console.log("MISSING figure: " + p); return; }
  children.push(new D.Paragraph({
    children: [new D.ImageRun({
      type: "png",
      data: fs.readFileSync(p),
      transformation: { width: 560, height: 400 },
    })],
    alignment: D.AlignmentType.CENTER,
    spacing: { before: 120, after: 80 },
  }));
  const legend = LEGEND[idx + 1] || ("**Fig. " + (idx + 1) + " | " + caption + "**");
  children.push(new D.Paragraph({
    children: [
      ...runs(legend),
      new D.TextRun({ text: "  [" + stem + ".pdf / .png]", italics: true, size: 18 }),
    ],
    spacing: { line: 240, after: 240 },
  }));
  if (idx < FIGS.length - 1)
    children.push(new D.Paragraph({ children: [new D.PageBreak()] }));
});
for (const t of FIGNOTE) {
  children.push(new D.Paragraph({ children: runs(t), spacing: { line: 240, after: 160 } }));
}

const doc = new D.Document({
  creator: "cqtna audit project",
  title: title || "Manuscript",
  styles: {
    default: {
      document: { run: { font: "Times New Roman", size: 24 } },   // 12pt
    },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal",
        quickFormat: true, run: { size: 28, bold: true, font: "Times New Roman" } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal",
        quickFormat: true, run: { size: 26, bold: true, font: "Times New Roman" } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal",
        quickFormat: true, run: { size: 24, bold: true, italics: true, font: "Times New Roman" } },
    ],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },        // US Letter
        margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 },
      },
      // lineNumbers belongs at the top of section properties, not inside page:
      // nested under page it is dropped silently and no w:lnNumType is written.
      lineNumbers: {
        countBy: 1,
        start: 1,
        restart: D.LineNumberRestartFormat.CONTINUOUS,
        distance: 360,
      },
    },
    headers: {
      default: new D.Header({
        children: [new D.Paragraph({
          alignment: D.AlignmentType.RIGHT,
          children: [new D.TextRun({ text: "Target nomination re-reads the outcome GWAS", size: 18, italics: true })],
        })],
      }),
    },
    footers: {
      default: new D.Footer({
        children: [new D.Paragraph({
          alignment: D.AlignmentType.CENTER,
          children: [new D.TextRun({ children: [D.PageNumber.CURRENT], size: 18 })],
        })],
      }),
    },
    children: [
      new D.Paragraph({
        children: runs(title || "Manuscript", { bold: true, size: 32 }),
        alignment: D.AlignmentType.CENTER,
        spacing: { after: 480 },
      }),
      ...children,
    ],
  }],
});

D.Packer.toBuffer(doc).then(buf => {
  const out = "D:/R_ex/MR/manuscript/MANUSCRIPT_GB.docx";
  fs.writeFileSync(out, buf);
  console.log("wrote " + out + " (" + (buf.length / 1024).toFixed(0) + " KB)");
});
