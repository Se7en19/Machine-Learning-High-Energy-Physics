# Simulación de muones y piones

## Configuración del detector. 

La simulación se basó en un módulo de detección compuesto por dos capas de material centellador plástico, las cuales estan organizadas de la siguiente manera: 

* **Composición:** Plástico centellador 
* **Estructura:** 1 módulo con 2 capas independientes.
  * **Capa 1:** 20 barras alieneadas y juntas
  * **Capa 2:** 20 barras perpendiculares a la primera capa 
* **Dimensiones de cada barra:** 
  * **Largo:** 1 metro 
  * **Ancho:** 5 centímetros. 
  * **Grosor:** 1 Centímetro (este parámetro es clave para las fluctuaciones de Landau en la pérdida de energía)

La disposición de las capas permite determinar la posición del impacto de la partícula en el plano "x, y". Por otro lado el grosor define el recorrido dentro de la partícula dentro del material sensible, lo cual determina la forma de pérdida de energía dE/dx observadas en los resultados.
Tambien es importante rresaltar que en la simulación se diaparan muones y piones tanto positivos como negativos, por lo que el resultado de estas simulaciones se basan en estas partículas. 

## Análisis de  los resultados
En esta sección, se presetan y explican las gráficas obtenidas tras procesar los datos  de la simulación en **ROOT**

### Idetificación por moemnto y carga ($p/q$)

![Separación PQ](../../Muon%20and%20Pion/img/dEdx_pq.png)

Al gráficar la pérdida de energía frente al momento dividido por la carga ($p/q$), podemos observar que, los dos brazos de la "V" representan las párticulas con carga positiva y negativa, lo que nos permite apreciar dos bandas distitntas, esto por otro lado, se debe a que el pion es más pesado que el muón, a un mismo momento tienen velocidades diferentes, lo que resulta en deósitos de energía distitntos. Gracias a esta separación, el detector es capaz de distinguir entre ambas partículas. 

### $dE/dx$ vs $\beta$

![dE/dx vs beta](../../Muon%20and%20Pion/img/dEdX_B.png)

En esta gráfica comparamos los resultados de la simulación con el modelo teórico de ***Bethe-Bloch** (la cual se representa por la línea negra discontinua). 
Podemos observar que a valores bajos de $\beta$ (partículas lentas), la pérdida de energía es altísima, A medida que se acercan a $\beta$ = 0.9, la pérdida de energía llega a un mínimo (el "MIP"). 
Esto quiere decir que en la simulación, al separar los paneles para $\mu$ y $\pi$, se miran sus nubes de pintos que son de forma casi identica, esto es físciamente correcto, porque $dE/dx$ dependen principalemnte de la velocidad y la carga.

### Dinámica de la pérdida de energía 

![dE/dx vs Momento](../../Muon%20and%20Pion/img/dEdx_momento.png)

Aquí se analiza la relación directa entre el momento ($p$) y la ionización en el material.
Se demuestra que, a medida que la partícula va más rápido, es decir hay momentos más altos (derecha),la curva se aplana, lo que quiere decir que las par´ticulas se comoprtan como "Partículas Mínimamente Ionizantes" (MIP), depositando la menor cantidad de energía posible, por otro laod, a momentos bajos (izquierda), la partícula pasa más tiempo cerca de los átomos del detector, transfiriendo más energía.

### $dE/dx$ vs $\beta$\gamma$

![dE/dx vs $\beta$\gamma$](../../Muon%20and%20Pion/img/dEdX_By.png)

ESta gráfica es la forma general de la pérdida de energía, ya que utiliza la variable $\beta$\gamma$. Al usar esta variable (el momento de la párticula dividido por la masa), las curvas de todas las partículas cargadas deberían colpasar en una sola forma. Ya que los datos de la gráfica siguen esta curva confirma que las propiedades el **plástico centellador** prueba que la simulación repordece interacciones reales. 

### Fluctuaciones de Landau Y energía cinética 

![Distribución de Landau](../../Muon%20and%20Pion/img/landau.png)

Esta distribución es el resultado de disparar partículas a través de las barras de 1 cm de grosor.
A diferencia de una campana de Gauss, la distribución de Landau tiene una "cola" hacia la derecha. Esto se debe a colisiones raras pero intensas que transfieren mucha energía a electrones individuales (rayos delta). Por lo que, el ancho de esta distribución esta ligado al grosor del plástico centellador definido en la geometría. 

### Identificación de partículas (ṔID) 

![PID](../../Muon%20and%20Pion/img/PID.png)

Esta gráfica presenta la comparación directa entre las diferentes especies de partículas simuladas.
Permite visualizar cómo se distribuyen los eventos de muones y piones en el espacio de fase de depósito de energía, esto nos ayuda a definir los cortes de selección (cuts) en el análisis de datos, determinando qué tan "limpia" puede ser la identificación de una partícula frente a otra en este diseño de detector centellador. 

### Análisis de Eficiencia en el sistema 

![Curva de Eficiencia](../../Muon%20and%20Pion/img/eficiencia.jpeg)

Finalmente, evaluamos el rendimiento del detector diseñado con barras de 1 metro de largo y 1 cm de grosor. La eficiencia se detemrina por: **número de eventos determinados /  número de eventos generados**. 
Se observa que el sistema mantiene una eficiencia cercana al **100%** para umbrales de energía por debajo de los 8 MeV, por otor lado la caída rápida después de los 10 MeV indica el límite físico del detector; las partículas raramente depositan más de esa energía en un trayecto de solo 1 cm, lo que ayuda a definir los parámetros de activación (trigger) del experimento real. 

