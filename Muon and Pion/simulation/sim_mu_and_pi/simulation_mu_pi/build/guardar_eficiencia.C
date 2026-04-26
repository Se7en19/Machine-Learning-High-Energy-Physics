void guardar_eficiencia() {
    // 1. Abrimos tu archivo en modo "UPDATE" (para poder meterle cosas nuevas)
    TFile *file = TFile::Open("output_run0.root", "UPDATE");
    if (!file || file->IsZombie()) {
        cout << "Error abriendo el archivo." << endl;
        return;
    }

    // 2. Leemos tus datos crudos
    TTree *tree = (TTree*)file->Get("Hits");
    Double_t fEdep;
    Int_t fEvent;
    tree->SetBranchAddress("fEdep", &fEdep);
    tree->SetBranchAddress("fEvent", &fEvent);

    // 3. Agrupamos la energía por evento (como lo hace un físico experimental)
    std::map<Int_t, Double_t> energia_por_evento;
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        energia_por_evento[fEvent] += fEdep;
    }

    int total_eventos = energia_por_evento.size();

    // 4. Calculamos la curva de eficiencia
    const int nPuntos = 100;
    Double_t umbrales[nPuntos];
    Double_t eficiencias[nPuntos];

    for (int i = 0; i < nPuntos; i++) {
        umbrales[i] = i * 2.0; // Revisamos umbrales de 0 a 200 MeV
        int detectados = 0;
        for (auto const& par : energia_por_evento) {
            if (par.second > umbrales[i]) detectados++;
        }
        eficiencias[i] = ((Double_t)detectados / total_eventos) * 100.0;
    }

    // 5. Creamos el objeto gráfico de ROOT
    TGraph *grafica = new TGraph(nPuntos, umbrales, eficiencias);
    grafica->SetName("Curva_Eficiencia"); // ¡ESTE ES EL NOMBRE QUE VERÁS EN EL TBROWSER!
    grafica->SetTitle("Eficiencia del Hodoscopio;Umbral (MeV);Eficiencia (%)");
    grafica->SetLineColor(kBlue);
    grafica->SetLineWidth(3);

    // 6. ¡LA MAGIA! Escribimos la gráfica dentro de tu output_run0.root
    grafica->Write();
    
    // Cerramos todo ordenadamente
    file->Close();
    cout << "¡Exito! La curva ha sido inyectada en tu archivo ROOT." << endl;
}
