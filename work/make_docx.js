// Gera o relatório .docx comparando os resultados 1.5 GL (QZS) x 2.5 GL (QZS-ADV).
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  ImageRun, Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
} = require("docx");

const IMG = "/home/user/Chaos/work";
const OUT = "/home/user/Chaos/Relatorio_QZS-ADV_1p5_vs_2p5GL.docx";

function img(file, wpx) {
  const dims = {
    "CMP_fracao_regular_preview.png": 2.75, "K5_absorvedor_preview.png": 1.5,
    "K3_FLI_preview.png": 1.341, "K7_FLI_2p5GL_preview.png": 1.098,
    "K2_secoes_preview.png": 1.25, "K8_secoes_2p5GL_preview.png": 1.25,
    "K6_amortecimento_preview.png": 1.786, "K9_amortecimento_2p5GL_preview.png": 1.25,
    "K4_fracao_preview.png": 1.4, "K1_esqueleto_preview.png": 1.833,
  };
  const asp = dims[file];
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 120, after: 60 },
    children: [new ImageRun({
      type: "png",
      data: fs.readFileSync(`${IMG}/${file}`),
      transformation: { width: wpx, height: Math.round(wpx / asp) },
    })],
  });
}

const cap = (t) => new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { after: 160 },
  children: [new TextRun({ text: t, italics: true, size: 18, color: "555555" })],
});
const P = (runs, opts = {}) => new Paragraph({ spacing: { after: 120 }, ...opts,
  children: Array.isArray(runs) ? runs : [new TextRun({ text: runs, size: 22 })] });
const T = (t, o = {}) => new TextRun({ text: t, size: 22, ...o });
const H1 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 240, after: 120 }, children: [new TextRun({ text: t, bold: true })] });
const H2 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 200, after: 100 }, children: [new TextRun({ text: t, bold: true })] });

// ---- tabelas ----
const COLS = [2600, 1160, 1900, 1900, 1800]; // soma 9360
function cell(text, { bold = false, fill = null, align = AlignmentType.LEFT, color = null } = {}, w) {
  return new TableCell({
    width: { size: w, type: WidthType.DXA },
    shading: fill ? { type: ShadingType.CLEAR, fill } : undefined,
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    children: [new Paragraph({ alignment: align, children: [new TextRun({ text, bold, size: 20, color: color || undefined })] })],
  });
}
function table(header, rows) {
  const hdr = new TableRow({
    tableHeader: true,
    children: header.map((h, i) => cell(h, { bold: true, fill: "1F4E8C", color: "FFFFFF", align: i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER }, COLS[i])),
  });
  const body = rows.map((r, ri) => new TableRow({
    children: r.map((v, i) => {
      const fill = ri % 2 ? "EEF2F8" : null;
      const isDelta = i === 4;
      const color = isDelta ? (v.startsWith("+") ? "1E7A1E" : "BC4B51") : null;
      return cell(v, { fill, align: i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER, color, bold: isDelta }, COLS[i]);
    }),
  }));
  return new Table({ columnWidths: COLS, width: { size: 9360, type: WidthType.DXA }, rows: [hdr, ...body] });
}

const regRows = [
  ["Monostable", "0,01", "100,0 %", "92,7 %", "-7,3"],
  ["Monostable", "0,05", "94,9 %", "45,9 %", "-49,0"],
  ["Monostable", "0,15", "39,0 %", "27,4 %", "-11,7"],
  ["Shallow wells", "0,01", "94,4 %", "82,3 %", "-12,1"],
  ["Shallow wells", "0,05", "79,3 %", "33,2 %", "-46,1"],
  ["Shallow wells", "0,15", "30,5 %", "19,6 %", "-10,8"],
  ["Deep wells", "0,01", "88,7 %", "78,2 %", "-10,6"],
  ["Deep wells", "0,05", "63,7 %", "61,5 %", "-2,2"],
  ["Deep wells", "0,15", "43,7 %", "40,9 %", "-2,9"],
];
const chaosRows = [
  ["Monostable", "0,002", "0,0 %", "12,5 %", "+12,5"],
  ["Monostable", "0,05", "15,6 %", "60,9 %", "+45,3"],
  ["Monostable", "0,15", "60,9 %", "76,6 %", "+15,6"],
  ["Shallow wells", "0,002", "9,4 %", "23,4 %", "+14,1"],
  ["Shallow wells", "0,05", "42,2 %", "73,4 %", "+31,2"],
  ["Shallow wells", "0,15", "81,2 %", "81,2 %", "+0,0"],
  ["Deep wells", "0,002", "4,7 %", "3,1 %", "-1,6"],
  ["Deep wells", "0,05", "20,3 %", "25,0 %", "+4,7"],
  ["Deep wells", "0,15", "34,4 %", "37,5 %", "+3,1"],
];

const children = [
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
    children: [new TextRun({ text: "Análise KAM do sistema QZS‑ADV", bold: true, size: 34 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
    children: [new TextRun({ text: "O que muda nos resultados ao incluir o absorvedor (1,5 GL → 2,5 GL)", size: 24, color: "444444" })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 240 },
    children: [new TextRun({ text: "Limite conservativo (ζ₁ = ζ₂ = 0) · μ = 0,1 · β = 0,35", italics: true, size: 20, color: "666666" })] }),

  H1("1. Contexto: os dois modelos"),
  P([
    T("Todas as figuras derivam do mesmo Hamiltoniano do oscilador QZS‑ADV. A análise foi feita em dois níveis. O "),
    T("modelo reduzido de 1,5 GL", { bold: true }),
    T(" congela o absorvedor (Z = W = 0): H₁ = p²/2 + U(X) − f·X·cos t, com U(X) = 1,55 X² − 3√((1,5η)² + X²). É o QZS puro. O "),
    T("modelo completo de 2,5 GL", { bold: true }),
    T(" (QZS‑ADV) integra o absorvedor acoplado: H = p²/2 + p₂²/(2μ) + U(X) + μβ²(X−X₂)²/2 − f·X·cos t, adicionando um grau de liberdade e um mapa tangente 4D."),
  ]),
  P([
    T("Antes, a análise de 2,5 GL existia apenas como um recorte (mapas FLI de "),
    T("shallow wells", { italics: true }),
    T("). Agora o sistema completo foi calculado para os "),
    T("três casos", { bold: true }),
    T(" (monostable, shallow wells, deep wells) em mapas FLI, seções de Poincaré e amortecimento. Este relatório resume o que muda nos resultados."),
  ]),

  H1("2. Resumo da mudança"),
  P([
    T("A conclusão central é direta: ", {}),
    T("acoplar o absorvedor reduz a fração regular (toros KAM) e aumenta o caos em praticamente todos os regimes", { bold: true }),
    T(". O grau de liberdade extra abre novos canais de ressonância que destroem toros que sobreviviam no modelo de 1,5 GL."),
  ]),
  img("CMP_fracao_regular_preview.png", 630),
  cap("Figura 1 — Fração regular do espaço de fase: QZS (1,5 GL, azul) vs QZS‑ADV (2,5 GL, vermelho), por caso e por amplitude de forçamento f."),
  P([
    T("Três padrões se destacam: (i) o efeito é ", {}),
    T("máximo em forçamento intermediário", { bold: true }),
    T(" (f = 0,05), onde a fração regular cai ~49 p.p. no monostable e ~46 p.p. no shallow wells; (ii) os poços profundos são os ", {}),
    T("mais robustos", { bold: true }),
    T(" (quedas de apenas 2–11 p.p.), por terem barreira alta e forte não linearidade; (iii) em f muito baixo o sistema ainda é quase todo regular nos dois modelos, e em f alto os dois já estão dominados pelo caos."),
  ]),

  H1("3. Números da comparação"),
  H2("3.1 Fração regular (mapas FLI, H₀ − U_min ≤ 0,5; FLI ≤ 8)"),
  table(["Caso", "f", "1,5 GL", "2,5 GL", "Δ (p.p.)"], regRows),
  P([T("A coluna Δ mostra a variação em pontos percentuais ao passar de 1,5 para 2,5 GL — sempre negativa: o absorvedor sempre encolhe a região regular.", { italics: true, size: 20, color: "555555" })], { spacing: { before: 80, after: 200 } }),
  H2("3.2 Fração caótica nas seções de Poincaré (FLI > 10)"),
  table(["Caso", "f", "1,5 GL", "2,5 GL", "Δ (p.p.)"], chaosRows),
  P([T("Espelho da tabela anterior: o caos cresce (Δ positivo) em quase todos os casos, com o maior salto no monostable a f = 0,05 (+45 p.p.). Os poços profundos mal se alteram.", { italics: true, size: 20, color: "555555" })], { spacing: { before: 80, after: 160 } }),

  H1("4. Comparação figura a figura"),

  H2("4.1 Mapas FLI — 1,5 GL vs 2,5 GL"),
  P("Nos mapas do indicador de Lyapunov rápido (escuro = toros KAM regulares; claro = caos), a mancha vermelha de caos cresce visivelmente ao ligar o absorvedor, sobretudo em f = 0,05."),
  img("K3_FLI_preview.png", 600),
  cap("Figura 2 — K3: mapas FLI do sistema QZS de 1,5 GL (absorvedor congelado)."),
  img("K7_FLI_2p5GL_preview.png", 470),
  cap("Figura 3 — K7: mapas FLI do sistema completo QZS‑ADV de 2,5 GL."),

  H2("4.2 Seções de Poincaré — 1,5 GL vs 2,5 GL"),
  P("As curvas invariantes coloridas (toros) aparecem mais fragmentadas e o mar de caos cinza é maior no modelo de 2,5 GL. Os poços profundos preservam bem os toros centrais em ambos."),
  img("K2_secoes_preview.png", 600),
  cap("Figura 4 — K2: seções de Poincaré do QZS de 1,5 GL."),
  img("K8_secoes_2p5GL_preview.png", 600),
  cap("Figura 5 — K8: seções de Poincaré do QZS‑ADV de 2,5 GL."),

  H2("4.3 Efeito direto do absorvedor"),
  P("A figura K5 sobrepõe os dois modelos para shallow wells: a linha de cima (1,5 GL) versus a de baixo (2,5 GL). A região regular encolhe de forma clara ao incluir o absorvedor — especialmente em f = 0,05 (79 % → 33 %)."),
  img("K5_absorvedor_preview.png", 520),
  cap("Figura 6 — K5: acoplamento ao absorvedor encolhe a região KAM (shallow wells)."),

  H2("4.4 Amortecimento — 1,5 GL vs 2,5 GL"),
  P([
    T("Fora do limite conservativo, o amortecimento destrói os toros e as órbitas colapsam em atratores. No modelo de 2,5 GL há dois coeficientes: ζ₁ (massa primária) e ζ₂ (absorvedor). ", {}),
    T("Observação de modelagem:", { bold: true }),
    T(" ζ₂ foi aplicado à velocidade absoluta do absorvedor; caso o modelo de referência use a velocidade relativa (típico de DVA), o cálculo é facilmente ajustável.", {}),
  ]),
  img("K6_amortecimento_preview.png", 600),
  cap("Figura 7 — K6: amortecimento no sistema de 1,5 GL (apenas ζ₁)."),
  img("K9_amortecimento_2p5GL_preview.png", 600),
  cap("Figura 8 — K9: amortecimento no sistema completo de 2,5 GL (pares ζ₁, ζ₂)."),

  H1("5. Conclusões"),
  P([T("• O absorvedor sempre reduz a regularidade — a inclusão do 2º grau de liberdade abre ressonâncias adicionais que quebram toros KAM.", {})]),
  P([T("• O efeito é máximo em forçamento intermediário (f ≈ 0,05) e nos regimes monostable e shallow wells; deep wells resistem.", {})]),
  P([T("• As seções de Poincaré e os mapas FLI concordam quantitativamente (mais caos ⇄ menos toros).", {})]),
  P([T("• A equivalência numérica com o cálculo original foi confirmada: a fração regular de shallow wells em 2,5 GL reproduz exatamente [0,823; 0,332; 0,196].", {})]),

  new Paragraph({ spacing: { before: 240 }, border: { top: { style: BorderStyle.SINGLE, size: 6, color: "CCCCCC", space: 8 } },
    children: [new TextRun({ text: "Nota: as imagens deste relatório são prévias geradas em Python (matplotlib) a partir de dados_kam.mat, para inspeção. As figuras finais de publicação (fundo transparente, PDF vetorial, 1900 dpi) são produzidas por fig_qzs_kam.m no MATLAB (fig_qzs_kam(1900)).", italics: true, size: 18, color: "777777" })] }),
];

const doc = new Document({
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1080, bottom: 1080, left: 1080, right: 1080 } } },
    children,
  }],
});
Packer.toBuffer(doc).then((b) => { fs.writeFileSync(OUT, b); console.log("escrito", OUT, b.length, "bytes"); });
