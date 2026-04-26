#include <iostream>
#include <set>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>

void analisis_eficiencia() {
    // =====================================================================
    // 1. CONFIGURACIÓN
    // =====================================================================
    double nGenerados = 10000.0; // Número de eventos por cada RUN
    int totalRuns = 10;          // Cambia esto al número de archivos que tengas
    
    // Aquí puedes poner los valores reales que cambiaste en tu barrido.
    // Ejemplo: si cambiaste la energía de 10 en 10 MeV: {10, 20, 30...}
    // Si no lo sabes, el eje X simplemente mostrará "0, 1, 2, 3..."
    double valoresX[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; 

    TGraph *grafica = new TGraph();
    grafica->SetTitle("Curva de Eficiencia del Detector;Variable de Barrido;Eficiencia (%)");

    // =====================================================================
    // 2. BUCLE PARA PROCESAR CADA ARCHIVO
    // =====================================================================
    for(int i = 0; i < totalRuns; i++) {
        TString nombreArchivo = Form("output_run%d.root", i);
        TFile *file = TFile::Open(nombreArchivo);
        
        if (!file || file->IsZombie()) {
            std::cout << "[-] Saltando: " << nombreArchivo << " (no encontrado)" << std::endl;
            continue;
        }

        TTree *tree = (TTree*)file->Get("Hits");
        if (!tree) {
            file->Close();
            continue;
        }

        int eventID;
        double edep;
        tree->SetBranchAddress("fEvent", &eventID);
        tree->SetBranchAddress("fEdep", &edep);

        // Usamos un set para contar eventos únicos por cada archivo
        std::set<int> eventosUnicos;

        for (long j = 0; j < tree->GetEntries(); j++) {
            tree->GetEntry(j);
            if (edep > 0.01) { // 0.01 MeV como umbral mínimo de detección
                eventosUnicos.insert(eventID);
            }
        }

        double eficiencia = ( (double)eventosUnicos.size() / nGenerados ) * 100.0;
        
        // Guardamos el punto en la gráfica
        // (i = eje X, eficiencia = eje Y)
        grafica->SetPoint(i, valoresX[i], eficiencia);

        std::cout << "[+] Archivo: " << nombreArchivo 
                  << " | Eficiencia: " << eficiencia << "%" << std::endl;

        file->Close();
    }

    // =====================================================================
    // 3. ESTILO Y DIBUJO
    // =====================================================================
    TCanvas *c1 = new TCanvas("c1", "Analisis de Eficiencia", 800, 600);
    c1->SetGrid();

    grafica->SetMarkerStyle(21);
    grafica->SetMarkerSize(1.2);
    grafica->SetMarkerColor(kBlue);
    grafica->SetLineColor(kBlue);
    grafica->SetLineWidth(2);

    grafica->Draw("ALP"); // A: ejes, L: línea, P: puntos
    
    // Forzamos el eje Y de 0 a 110 para ver bien la curva
    grafica->GetYaxis()->SetRangeUser(0, 110);

    c1->SaveAs("grafica_eficiencia_completa.png");
}
