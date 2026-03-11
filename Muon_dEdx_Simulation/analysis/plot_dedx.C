// ============================================================================
// plot_dedx.C
// Macro ROOT para análisis de la simulación de muones en Geant4 B2a
//
// Genera 6 gráficas:
//   --- Bethe-Bloch ---
//   1. dE/dx vs p/q        (rigidez magnética — separa μ⁺ y μ⁻)
//   2. dE/dx vs β          (velocidad reducida)
//   3. dE/dx vs βγ         (variable estándar en papers — eje log)
//   --- Hits ---
//   4. Hits por cámara     (distribución de chamberNb)
//   5. Distribución Landau (edep por cámara — una por cada una de las 5)
//   6. Mapa espacial       (posición x vs y de los hits)
//
// Ejecutar desde build/:
//   root -l '../analysis/plot_dedx.C'
// ============================================================================

// ── Parámetros físicos del Xenón gaseoso ────────────────────────────────────
namespace XeParams {
    const double Z   = 54.0;
    const double A   = 131.29;
    const double I   = 482.0e-6;   // potencial de ionización [MeV]
    const double rho = 2.9e-3;     // densidad efectiva [g/cm³]
    const double K   = 0.307075;   // constante de Bethe-Bloch [MeV cm²/mol]
    const double me  = 0.511;      // masa electrón [MeV/c²]
    const double mmu = 105.6584;   // masa muón [MeV/c²]
}

// ── Bethe-Bloch: recibe βγ, devuelve dE/dx en MeV/mm ───────────────────────
double BetheBloch(double bg)
{
    using namespace XeParams;
    if (bg <= 0.) return 0.;

    double beta2 = bg*bg / (1.0 + bg*bg);
    double gamma = std::sqrt(1.0 + bg*bg);

    double Tmax = (2.0 * me * bg*bg) /
                  (1.0 + 2.0*gamma*(me/mmu) + (me/mmu)*(me/mmu));
    if (Tmax <= 0.) return 0.;

    double logArg = (2.0 * me * beta2 * gamma*gamma * Tmax) / (I * I);
    if (logArg <= 1.) return 0.;

    double dedx_cm = K * (Z/A) * (1.0/beta2) *
                     (0.5 * std::log(logArg) - beta2);
    return dedx_cm * rho / 10.0;
}

// ── Estética global ─────────────────────────────────────────────────────────
void SetStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetTickLength(0.02, "XY");
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.16);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadTopMargin(0.09);
}

// ── Canvas con resolución fija alta (1600×1200 px) ──────────────────────────
// Equivale a ~200 DPI en una figura de 8×6 pulgadas (tamaño estándar de paper)
// ROOT guarda PNG a 72 DPI por defecto, pero al crear el canvas más grande
// el PNG resultante tiene más píxeles y se ve nítido en cualquier pantalla.
TCanvas* MakeCanvas(const char* name, const char* title,
                    bool leftExtra = false)
{
    // 1600 × 1200 píxeles — resolución fija para todas las gráficas
    TCanvas* c = new TCanvas(name, title, 1600, 1200);

    // Márgenes estándar
    c->SetLeftMargin(leftExtra ? 0.16 : 0.13);
    c->SetRightMargin(0.16);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.09);
    return c;
}

// ── Etiqueta informativa superior ───────────────────────────────────────────
void AddInfoLabel(TCanvas*, const char* extra = "")
{
    TLatex* info = new TLatex();
    info->SetNDC();
    info->SetTextFont(42);
    info->SetTextSize(0.032);  // ligeramente más grande para 1600px
    info->SetTextColor(kBlack);
    TString label = TString::Format(
        "Geant4 B2a  |  #mu^{#pm} in Xe (gas)  |  "
        "FTFP\\_BERT  |  10^{5} events  %s", extra);
    info->DrawLatex(0.13, 0.938, label.Data());
}

// ── Estética de histograma 2D ────────────────────────────────────────────────
void StyleHisto2D(TH2D* h, const char* xt, const char* yt)
{
    h->GetXaxis()->SetTitle(xt);
    h->GetYaxis()->SetTitle(yt);
    h->GetZaxis()->SetTitle("Counts");
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetZaxis()->SetTitleSize(0.038);
    h->GetXaxis()->SetLabelSize(0.038);
    h->GetYaxis()->SetLabelSize(0.038);
    h->GetZaxis()->SetLabelSize(0.034);
    h->GetXaxis()->SetTitleOffset(1.10);
    h->GetYaxis()->SetTitleOffset(1.50);
    h->GetZaxis()->SetTitleOffset(1.30);
}

// ── Estética de histograma 1D ────────────────────────────────────────────────
void StyleHisto1D(TH1* h, const char* xt, const char* yt,
                  Color_t fill = kAzure-3, Color_t line = kBlue+2)
{
    h->GetXaxis()->SetTitle(xt);
    h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.038);
    h->GetYaxis()->SetLabelSize(0.038);
    h->GetXaxis()->SetTitleOffset(1.10);
    h->GetYaxis()->SetTitleOffset(1.50);
    h->SetFillColorAlpha(fill, 0.6);
    h->SetLineColor(line);
    h->SetLineWidth(2);
}

// ============================================================================
// GRÁFICA 1: dE/dx vs p/q
// ============================================================================
void Plot1_pq(TTree* tree)
{
    TH2D* h = new TH2D("h_pq", "",
                        400, -3.2,    3.2,
                        200,  0.0003, 0.0010);

    tree->Draw("dedx:(rigidity/1000.0)>>h_pq",
               "dedx>0.0003 && dedx<0.0010 && "
               "TMath::Abs(rigidity/1000.0)>0.1 && "
               "TMath::Abs(rigidity/1000.0)<3.2",
               "COLZ goff");

    TCanvas* c = MakeCanvas("c_pq", "dE/dx vs p/q", true);
    c->SetLogz();

    StyleHisto2D(h, "p/q  (GeV/c)", "dE/dx  (MeV/mm)");
    h->GetYaxis()->SetTitleOffset(1.40);
    h->Draw("COLZ");

    // Bethe-Bloch lado positivo
    TF1* bb_pos = new TF1("bb_pos",
        [](double* x, double*) -> double {
            if (x[0] < 0.1) return 0.;
            return BetheBloch(TMath::Abs(x[0])*1000.0 / XeParams::mmu);
        }, 0.1, 3.2, 0);
    bb_pos->SetNpx(2000);
    bb_pos->SetLineColor(kBlack);
    bb_pos->SetLineWidth(2);
    bb_pos->Draw("same");

    // Bethe-Bloch lado negativo
    TF1* bb_neg = new TF1("bb_neg",
        [](double* x, double*) -> double {
            if (x[0] > -0.1) return 0.;
            return BetheBloch(TMath::Abs(x[0])*1000.0 / XeParams::mmu);
        }, -3.2, -0.1, 0);
    bb_neg->SetNpx(2000);
    bb_neg->SetLineColor(kBlack);
    bb_neg->SetLineWidth(2);
    bb_neg->Draw("same");

    // Línea vertical en 0
    TLine* line = new TLine(0., 0.0003, 0., 0.0010);
    line->SetLineColor(kGray+1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    // Etiquetas μ⁺ μ⁻
    TLatex* lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(42);
    lat->SetTextSize(0.045);
    lat->SetTextColor(kWhite);
    lat->DrawLatex(0.18, 0.75, "#mu^{-}");
    lat->DrawLatex(0.68, 0.75, "#mu^{+}");

    AddInfoLabel(c);
    c->Update();
    c->SaveAs("dedx_vs_pq.png");
    c->SaveAs("dedx_vs_pq.pdf");
    std::cout << "[1] Saved: dedx_vs_pq.png/.pdf  ("
              << h->GetEntries() << " entries)" << std::endl;
}

// ============================================================================
// GRÁFICA 2: dE/dx vs β
// ============================================================================
void Plot2_beta(TTree* tree)
{
    const char* beta_expr =
        "TMath::Abs(rigidity)/TMath::Sqrt(rigidity*rigidity+105.6584*105.6584)";

    TH2D* h = new TH2D("h_beta", "",
                        300, 0.3,    1.0,
                        200, 0.0003, 0.0010);

    TString draw = TString::Format("dedx:(%s)>>h_beta", beta_expr);
    tree->Draw(draw.Data(),
               "dedx>0.0003 && dedx<0.0010 && "
               "TMath::Abs(rigidity)>100 && TMath::Abs(rigidity)<3200",
               "COLZ goff");

    TCanvas* c = MakeCanvas("c_beta", "dE/dx vs beta", true);
    c->SetLogz();

    StyleHisto2D(h, "#beta = v/c", "dE/dx  (MeV/mm)");
    h->GetYaxis()->SetTitleOffset(1.40);
    h->Draw("COLZ");

    TF1* bb = new TF1("bb_beta",
        [](double* x, double*) -> double {
            double beta = x[0];
            if (beta <= 0. || beta >= 1.) return 0.;
            double bg = beta / std::sqrt(1.0 - beta*beta);
            return BetheBloch(bg);
        }, 0.3, 1.0, 0);
    bb->SetNpx(2000);
    bb->SetLineColor(kBlack);
    bb->SetLineWidth(2);
    bb->Draw("same");

    // Línea MIP en βγ ≈ 3.5 → β ≈ 0.962
    double beta_mip = 3.5 / std::sqrt(1.0 + 3.5*3.5);
    TLine* mip = new TLine(beta_mip, 0.0003, beta_mip, 0.0010);
    mip->SetLineColor(kYellow+1);
    mip->SetLineStyle(2);
    mip->SetLineWidth(2);
    mip->Draw("same");

    TLatex* lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(42);
    lat->SetTextSize(0.038);
    lat->SetTextColor(kWhite);
    lat->DrawLatex(0.20, 0.75, "#mu^{+} + #mu^{-}  (superimposed)");

    TLatex* mip_lat = new TLatex();
    mip_lat->SetNDC();
    mip_lat->SetTextFont(42);
    mip_lat->SetTextSize(0.030);
    mip_lat->SetTextColor(kYellow+1);
    mip_lat->DrawLatex(0.58, 0.18, "MIP  (#beta#gamma #approx 3.5)");

    AddInfoLabel(c);
    c->Update();
    c->SaveAs("dedx_vs_beta.png");
    c->SaveAs("dedx_vs_beta.pdf");
    std::cout << "[2] Saved: dedx_vs_beta.png/.pdf  ("
              << h->GetEntries() << " entries)" << std::endl;
}

// ============================================================================
// GRÁFICA 3: dE/dx vs βγ
// ============================================================================
void Plot3_betagamma(TTree* tree)
{
    const char* bg_expr = "TMath::Abs(rigidity)/105.6584";

    TH2D* h = new TH2D("h_bg", "",
                        300, 0.5,    30.0,
                        200, 0.0003, 0.0010);

    TString draw = TString::Format("dedx:(%s)>>h_bg", bg_expr);
    tree->Draw(draw.Data(),
               "dedx>0.0003 && dedx<0.0010 && "
               "TMath::Abs(rigidity)>50 && TMath::Abs(rigidity)<3200",
               "COLZ goff");

    TCanvas* c = MakeCanvas("c_bg", "dE/dx vs betagamma", true);
    c->SetLogz();
    c->SetLogx();

    StyleHisto2D(h, "#beta#gamma = p/mc", "dE/dx  (MeV/mm)");
    h->GetXaxis()->SetMoreLogLabels();
    h->Draw("COLZ");

    TF1* bb = new TF1("bb_bg",
        [](double* x, double*) -> double {
            return BetheBloch(x[0]);
        }, 0.5, 30.0, 0);
    bb->SetNpx(2000);
    bb->SetLineColor(kBlack);
    bb->SetLineWidth(2);
    bb->Draw("same");

    // Línea MIP
    TLine* mip = new TLine(3.5, 0.0003, 3.5, 0.0010);
    mip->SetLineColor(kYellow+1);
    mip->SetLineStyle(2);
    mip->SetLineWidth(2);
    mip->Draw("same");

    TLatex lat;
    lat.SetTextFont(42);
    lat.SetTextSize(0.032);
    lat.SetTextColor(kWhite);
    lat.DrawLatex(0.6, 0.00090, "1/#beta^{2} region");
    lat.DrawLatex(15.0, 0.00058, "Rel. rise");

    TLatex mip_lat;
    mip_lat.SetNDC();
    mip_lat.SetTextFont(42);
    mip_lat.SetTextSize(0.028);
    mip_lat.SetTextColor(kYellow+1);
    mip_lat.DrawLatex(0.50, 0.18, "MIP");

    TLatex mu_lat;
    mu_lat.SetNDC();
    mu_lat.SetTextFont(42);
    mu_lat.SetTextSize(0.040);
    mu_lat.SetTextColor(kWhite);
    mu_lat.DrawLatex(0.20, 0.75, "#mu^{+} + #mu^{-}");

    AddInfoLabel(c);
    c->Update();
    c->SaveAs("dedx_vs_betagamma.png");
    c->SaveAs("dedx_vs_betagamma.pdf");
    std::cout << "[3] Saved: dedx_vs_betagamma.png/.pdf  ("
              << h->GetEntries() << " entries)" << std::endl;
}

// ============================================================================
// GRÁFICA 4: Hits por cámara
// Muestra cuántos hits registró cada una de las 5 cámaras.
// Útil para verificar que la geometría es correcta y todas las
// cámaras son igualmente eficientes.
// ============================================================================
void Plot4_hits_per_chamber(TTree* tree)
{
    // 5 bins, uno por cámara (0-4)
    TH1D* h = new TH1D("h_chambers", "", 5, -0.5, 4.5);

    tree->Draw("chamberNb>>h_chambers", "", "goff");

    TCanvas* c = MakeCanvas("c_chambers", "Hits per chamber", false);

    StyleHisto1D(h, "Chamber number", "Number of hits");

    // Forzar eje Y desde 0 para que las diferencias entre cámaras sean reales
    h->SetMinimum(0);
    double ymax = h->GetMaximum();
    h->SetMaximum(ymax * 1.20);  // 20% de espacio arriba para las etiquetas

    // Etiquetas en el eje X: 0, 1, 2, 3, 4
    for (int i = 0; i < 5; i++) {
        h->GetXaxis()->SetBinLabel(i+1, TString::Format("%d", i).Data());
    }
    h->GetXaxis()->SetLabelSize(0.050);

    h->Draw("HIST");

    // Valor numérico sobre cada barra
    for (int i = 1; i <= 5; i++) {
        TLatex* val = new TLatex();
        val->SetTextFont(42);
        val->SetTextSize(0.032);
        val->SetTextAlign(22);
        val->DrawLatex(h->GetBinCenter(i),
                       h->GetBinContent(i) * 1.03,
                       TString::Format("%.0f", h->GetBinContent(i)).Data());
    }

    // Línea de media para referencia
    double mean = h->GetMean(2);  // media en Y
    TLine* avg = new TLine(-0.5, mean, 4.5, mean);
    avg->SetLineColor(kRed);
    avg->SetLineStyle(2);
    avg->SetLineWidth(2);
    avg->Draw("same");

    TLatex* avg_lat = new TLatex();
    avg_lat->SetNDC();
    avg_lat->SetTextFont(42);
    avg_lat->SetTextSize(0.028);
    avg_lat->SetTextColor(kRed);
    avg_lat->DrawLatex(0.70, 0.55,
        TString::Format("Mean = %.0f hits/chamber", mean).Data());

    // Nota explicativa
    TLatex* note = new TLatex();
    note->SetNDC();
    note->SetTextFont(42);
    note->SetTextSize(0.026);
    note->SetTextColor(kGray+2);
    note->DrawLatex(0.13, 0.13,
        "Each hit = one step with E_{dep} > 0 inside the Xe chamber");

    AddInfoLabel(c);
    c->Update();
    c->SaveAs("hits_per_chamber.png");
    c->SaveAs("hits_per_chamber.pdf");
    std::cout << "[4] Saved: hits_per_chamber.png/.pdf" << std::endl;

    // Imprimir estadísticas en consola
    std::cout << "    Chamber hits: ";
    for (int i = 1; i <= 5; i++)
        std::cout << "Ch" << i-1 << "=" << (int)h->GetBinContent(i) << "  ";
    std::cout << std::endl;
}

// ============================================================================
// GRÁFICA 5: Distribución de Landau — edep por cámara
// Muestra la distribución de energía depositada en cada paso dentro
// de cada cámara. La forma asimétrica con cola larga es la distribución
// de Landau, característica de la ionización en capas delgadas.
// ============================================================================
void Plot5_landau(TTree* tree)
{
    // Colores para las 5 cámaras
    const Color_t colors[5] = {
        kBlue+1, kRed+1, kGreen+2, kOrange+1, kViolet+1
    };
    const Color_t fills[5] = {
        kBlue-9, kRed-9, kGreen-9, kOrange-9, kViolet-9
    };

    // Encontrar el rango real de edep primero
    double edep_max = tree->GetMaximum("dedx") * 0.8;  // corte en 80% del max

    TCanvas* c = MakeCanvas("c_landau", "Landau distribution per chamber", false);
    c->SetLogy();  // eje Y logarítmico para ver la cola de Landau

    TLegend* leg = new TLegend(0.65, 0.55, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.032);

    TH1D* hists[5];
    double global_max = 0;

    // Crear y llenar un histograma por cámara
    for (int ch = 0; ch < 5; ch++) {
        hists[ch] = new TH1D(
            TString::Format("h_landau_ch%d", ch), "",
            150, 0., 0.008);  // rango en MeV/mm

        TString draw = TString::Format("dedx>>h_landau_ch%d", ch);
        TString cut  = TString::Format(
            "chamberNb==%d && dedx>0 && dedx<0.008", ch);
        tree->Draw(draw.Data(), cut.Data(), "goff");

        hists[ch]->SetLineColor(colors[ch]);
        hists[ch]->SetLineWidth(2);
        hists[ch]->SetFillColorAlpha(fills[ch], 0.35);

        if (hists[ch]->GetMaximum() > global_max)
            global_max = hists[ch]->GetMaximum();

        leg->AddEntry(hists[ch],
            TString::Format("Chamber %d  (MPV = %.4f MeV/mm)",
                ch, hists[ch]->GetXaxis()->GetBinCenter(
                    hists[ch]->GetMaximumBin())).Data(),
            "lf");
    }

    // Dibujar — primero el que tiene mayor amplitud para fijar el rango
    hists[0]->GetXaxis()->SetTitle("dE/dx per step  (MeV/mm)");
    hists[0]->GetYaxis()->SetTitle("Entries");
    hists[0]->GetXaxis()->SetTitleSize(0.050);
    hists[0]->GetYaxis()->SetTitleSize(0.050);
    hists[0]->GetXaxis()->SetLabelSize(0.040);
    hists[0]->GetYaxis()->SetLabelSize(0.040);
    hists[0]->GetYaxis()->SetTitleOffset(1.10);
    hists[0]->SetMaximum(global_max * 5.0);  // espacio para eje log
    hists[0]->SetMinimum(0.5);
    hists[0]->Draw("HIST");

    for (int ch = 1; ch < 5; ch++)
        hists[ch]->Draw("HIST same");

    leg->Draw();

    // Nota sobre la distribución de Landau
    TLatex* note = new TLatex();
    note->SetNDC();
    note->SetTextFont(42);
    note->SetTextSize(0.026);
    note->SetTextColor(kGray+2);
    note->DrawLatex(0.13, 0.13,
        "Landau distribution: asymmetric tail from high-energy #delta-rays");

    // Flecha apuntando a la cola
    TArrow* arr = new TArrow(0.005, global_max*0.5,
                              0.006, global_max*0.08,
                              0.012, ">");
    arr->SetLineColor(kGray+2);
    arr->SetFillColor(kGray+2);
    arr->Draw();

    TLatex* tail = new TLatex(0.0045, global_max*0.6, "Landau tail");
    tail->SetTextFont(42);
    tail->SetTextSize(0.028);
    tail->SetTextColor(kGray+2);
    tail->Draw();

    AddInfoLabel(c);
    c->Update();
    c->SaveAs("landau_per_chamber.png");
    c->SaveAs("landau_per_chamber.pdf");
    std::cout << "[5] Saved: landau_per_chamber.png/.pdf" << std::endl;
}

// ============================================================================
// GRÁFICA 6: Mapa espacial de hits — x vs y
// Muestra dónde ocurrieron los hits en el plano transverso (x-y).
// Debe verse la sección circular de las cámaras de Xenón.
// Los hits del blanco de Plomo NO aparecen (no es volumen sensible).
// ============================================================================
void Plot6_spatial(TTree* tree)
{
    // Verificar si la rama pos existe en el TTree
    // (requiere que TrackerHit guarde la posición x,y,z)
    // Si no existe, usa chamberNb como proxy de posición z

    TBranch* bx = tree->GetBranch("x");
    TBranch* by = tree->GetBranch("y");

    TCanvas* c = MakeCanvas("c_spatial", "Spatial hit map", true);

    if (bx && by) {
        // ── Caso A: tenemos coordenadas x, y almacenadas ────────────────
        TH2D* h = new TH2D("h_xy", "",
                            200, -1500., 1500.,   // x en mm
                            200, -1500., 1500.);  // y en mm

        tree->Draw("y:x>>h_xy", "dedx>0", "COLZ goff");

        gStyle->SetPalette(kRainBow);
        StyleHisto2D(h, "x  (mm)", "y  (mm)");
        h->Draw("COLZ");

        // Círculo indicando el radio máximo del tracker
        TEllipse* tracker = new TEllipse(0., 0., 1200., 1200.);
        tracker->SetLineColor(kWhite);
        tracker->SetLineStyle(2);
        tracker->SetLineWidth(2);
        tracker->SetFillStyle(0);
        tracker->Draw("same");

        TLatex* rad_lat = new TLatex();
        rad_lat->SetNDC();
        rad_lat->SetTextFont(42);
        rad_lat->SetTextSize(0.028);
        rad_lat->SetTextColor(kWhite);
        rad_lat->DrawLatex(0.55, 0.87, "Tracker boundary");

    } else {
        // ── Caso B: no tenemos x,y → graficar chamberNb vs dedx ────────
        c->SetLogz();
        c->SetLeftMargin(0.16);
        c->SetRightMargin(0.16);

        // Rango Y ajustado al rango real de los datos (sin los outliers extremos)
        // Usar solo dedx < 0.003 para mostrar la región principal con claridad
        TH2D* h = new TH2D("h_ch_dedx", "",
                            5,   -0.5,   4.5,
                            200,  1e-4,  3e-3);  // rango logarítmico efectivo

        tree->Draw("dedx:chamberNb>>h_ch_dedx",
                   "dedx>1e-4 && dedx<3e-3", "COLZ goff");

        StyleHisto2D(h, "Chamber number", "dE/dx per step  (MeV/mm)");
        h->GetYaxis()->SetTitleOffset(1.40);

        // Etiquetas eje X
        for (int i = 0; i < 5; i++)
            h->GetXaxis()->SetBinLabel(i+1,
                TString::Format("Ch %d", i).Data());
        h->GetXaxis()->SetLabelSize(0.045);

        h->Draw("COLZ");

        // Línea de MPV
        TLine* mpv = new TLine(-0.5, 4e-4, 4.5, 4e-4);
        mpv->SetLineColor(kWhite);
        mpv->SetLineStyle(2);
        mpv->SetLineWidth(2);
        mpv->Draw("same");

        TLatex* mpv_lat = new TLatex();
        mpv_lat->SetNDC();
        mpv_lat->SetTextFont(42);
        mpv_lat->SetTextSize(0.030);
        mpv_lat->SetTextColor(kWhite);
        mpv_lat->DrawLatex(0.15, 0.32, "Landau MPV (#approx 4#times10^{-4} MeV/mm)");

        TLatex* note = new TLatex();
        note->SetNDC();
        note->SetTextFont(42);
        note->SetTextSize(0.026);
        note->SetTextColor(kGray+2);
        note->DrawLatex(0.13, 0.13,
            "x,y branches not stored — showing dE/dx distribution per chamber");
    }

    AddInfoLabel(c);
    c->Update();
    c->SaveAs("spatial_hit_map.png");
    c->SaveAs("spatial_hit_map.pdf");
    std::cout << "[6] Saved: spatial_hit_map.png/.pdf" << std::endl;
}

// ============================================================================
// FUNCIÓN PRINCIPAL
// ============================================================================
void plot_dedx()
{
    // ── Abrir archivo ──────────────────────────────────────────────────
    TFile* f = TFile::Open("B2a_output.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open B2a_output.root" << std::endl;
        return;
    }

    TTree* tree = (TTree*)f->Get("Hits");
    if (!tree) {
        std::cerr << "Error: TTree 'Hits' not found" << std::endl;
        f->ls();
        return;
    }

    std::cout << "=====================================================" << std::endl;
    std::cout << " Muon dE/dx Analysis — Geant4 B2a"                    << std::endl;
    std::cout << " Total entries: " << tree->GetEntries()               << std::endl;
    std::cout << "=====================================================" << std::endl;

    // ── Mostrar ramas disponibles ──────────────────────────────────────
    std::cout << "\nAvailable branches:" << std::endl;
    tree->Print();

    // ── Estética global ────────────────────────────────────────────────
    SetStyle();

    // ── Bethe-Bloch plots ──────────────────────────────────────────────
    std::cout << "\n--- Bethe-Bloch plots ---" << std::endl;
    Plot1_pq(tree);
    Plot2_beta(tree);
    Plot3_betagamma(tree);

    // ── Hit plots ──────────────────────────────────────────────────────
    std::cout << "\n--- Hit analysis plots ---" << std::endl;
    Plot4_hits_per_chamber(tree);
    Plot5_landau(tree);
    Plot6_spatial(tree);

    // ── Resumen final ──────────────────────────────────────────────────
    std::cout << "\n=====================================================" << std::endl;
    std::cout << " All 6 plots generated:"                               << std::endl;
    std::cout << "  [Bethe-Bloch]"                                       << std::endl;
    std::cout << "    dedx_vs_pq.png/.pdf"                               << std::endl;
    std::cout << "    dedx_vs_beta.png/.pdf"                             << std::endl;
    std::cout << "    dedx_vs_betagamma.png/.pdf"                        << std::endl;
    std::cout << "  [Hits]"                                              << std::endl;
    std::cout << "    hits_per_chamber.png/.pdf"                         << std::endl;
    std::cout << "    landau_per_chamber.png/.pdf"                       << std::endl;
    std::cout << "    spatial_hit_map.png/.pdf"                          << std::endl;
    std::cout << "=====================================================" << std::endl;
}
