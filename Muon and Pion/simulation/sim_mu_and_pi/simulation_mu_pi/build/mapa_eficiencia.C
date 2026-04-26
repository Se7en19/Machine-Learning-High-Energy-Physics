void mapa_eficiencia() {
    // 1. Abrimos tu archivo
    TFile *file = TFile::Open("output_run0.root");
    if (!file || file->IsZombie()) return;

    TTree *tree = (TTree*)file->Get("Hits");
    
    Double_t fX, fY, fEdep;
    tree->SetBranchAddress("fX", &fX);
    tree->SetBranchAddress("fY", &fY);
    tree->SetBranchAddress("fEdep", &fEdep);

    // 2. Creamos dos mapas vacios (de -500 a 500 mm)
    // h_total: Cuenta TODOS los impactos
    // h_detectado: Cuenta solo los que superan el umbral
    TH2F *h_total = new TH2F("h_total", "Total", 50, -500, 500, 50, -500, 500);
    TH2F *h_detectado = new TH2F("h_detectado", "Detectados", 50, -500, 500, 50, -500, 500);

    // Definimos tu umbral de energia (Ej. 1 MeV)
    Double_t umbral = 1.0; 

    // 3. Llenamos los mapas
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        h_total->Fill(fX, fY); // Aqui contamos todos
        
        if (fEdep > umbral) {
            h_detectado->Fill(fX, fY); // Aqui solo los exitosos
        }
    }

    // 4. ¡Matematicas! Dividimos los exitosos entre el total para sacar la Eficiencia
    TH2F *h_eficiencia = (TH2F*)h_detectado->Clone("h_eficiencia");
    h_eficiencia->Divide(h_total);
    h_eficiencia->Scale(100.0); // Lo multiplicamos por 100 para que sea Porcentaje (%)

    // 5. A dibujar con colores
    TCanvas *c1 = new TCanvas("c1", "Mapa de Eficiencia", 800, 600);
    c1->SetRightMargin(0.15); // Espacio para la barra de colores
    gStyle->SetPalette(kRainBow); // Paleta arcoiris termica
    gStyle->SetOptStat(0); // Escondemos el cuadro feo de estadisticas

    h_eficiencia->SetTitle("Mapa de Eficiencia Espacial del Hodoscopio;Posicion X (mm);Posicion Y (mm);Eficiencia (%)");
    
    // colz dibuja la matriz de colores y la barra lateral
    h_eficiencia->Draw("colz"); 

    // Guardamos tu obra de arte
    c1->SaveAs("mi_eficiencia_termica.png");
}
