# multimodal_posturography
Description of analysis of kinetic data of center of pressure (COP) from a ASCII file of force platform assessment of sensory organization test.

#on SciLab

#open file
#read header
    #name iniciais, year of birthday, ser, date of test, test, conditions, trials 
    #first tab 
#kinetic measures
    #idetify test, frequency of execution, data available, select variables defined
    #odentify calibration if avaiable
    #identify load cell forces (5)
    #confirm force plate components (6?)
    #calculate Fx, Fy, Fz (x - ap, y - ml, z - vertical for Marcos Duarte; on NeuroCom x - ml, y - ap, ver z)
    #calculate Mx, My, Mz
    #calculate CP (position defined by two coordinates on the platform surface according to the orientation of the individual) - calculados como CPap=(-h*Fx −My)/Fz e CPml=(−h*Fy+Mx)/Fz, em que h é a altura da base de apoio acima da plataforma de força por exemplo, um tapete sobre a plataforma de força.
    # Dados do CoP: estatocinesigrama e estabilograma. O estatocinesigrama é o mapa do CP na direção ap versus o CP na direção ml, o estabilograma é a série temporal do CP em cada uma das direções: ap e ml.
    // A posição do CG é uma medida de deslocamento e é totalmente independente da velocidade ou aceleração total do corpo e de seus segmentos. O CP também é uma medida de deslocamento e é dependente do CG, mas o CP expressa a localização do vetor resultante da força de reação do solo em uma plataforma de força. Esse vetor é igual e oposto à média ponderada da localização de todas as forças que agem na plataforma de força, como a força peso e as forças internas (musculares e articulares) transmitidas ao chão. O deslocamento do CG é a grandeza que realmente indica a oscilação do corpo inteiro, e a grandeza CP, é na verdade, uma combinação da resposta neuromuscular ao deslocamento do CG e da própria posição do CG.
    // As diferenças entre o CG e o CP são relacionadas à aceleração do corpo e, quanto menores as frequências de oscilação do corpo,  menores serão as diferenças entre essas duas grandezas.
    // O componente do CG numa direção horizontal é denominado projeção vertical do CG (CGv).
    #calculate projeção vertical do CG (CGv)
    #calculate a diferença de CP e CGv (CP-CGv)
    // o componente horizontal do CG, o CGv, pode ser estimada por integração dupla da força horizontal dividida pela massa (aceleração horizontal). Achar a posição e a velocidade inicial do corpo após a dupla integração. constantes. Considerar que no instante em que a força horizontal é nula, as posições do CP e da CGv coincidem. Zatsiorsky e Duarte propuseram considerando ambas as constantes de integração determinadas analiticamente a partir dos dados do CP, e os instantes de força nula são determinados pela interpolação dos dados obtidos da série temporal do CP.
    // Pêndulo invertido: A estimativa da CGv acontece a partir do CP, usando o método de filtragem a partir da relação, no domínio de frequências, entre CP e CGv, considerando o corpo como um pêndulo invertido. Esse método consiste na aplicação de um filtro passa-baixa no sinal do CP. A frequência de corte desse filtro é determinada a partir das características antropométricas do corpo, e geralmente a frequência é da ordem de 0,5 Hz.
    // estacionariedade do CoP - filtro passa-alta
    #filtragem - filtro passa-baixa de 10 Hz ou mais. A frequência do filtro deve ser escolhida em função de parâmetros da tarefa e do equipamento utilizado.
    //As variáveis podem ser derivadas do estatocinesigrama e estabilograma do CP. Classes da análise: análise global e análise estrutural. A análise global está relacionada à mensuração do “tamanho” dos padrões de oscilação tanto no domínio do tempo como no domínio das frequências. A análise estrutural identifica subunidades nos dados posturográficos e as relaciona aos processos de controle motor. Análise global: a trajetória do CP e a banda de frequência do estabilograma; análise estrutural: a velocidade média do
CP tem sido considerada a medida com maior confiabilidade entre as repetições.
// Códigos partem do pressuposto de que os dados do CP nas direções ap e ml, respectivamente como CPap e CPml, são variáveis no ambiente  Matlab.
//Análise global: a posição média do CP não é de interesse, pois ela é simplesmente dependente da posição absoluta do sujeito sobre a plataforma de força, a qual, geralmente, não é controlada. Nesse sentido, é um procedimento comum remover a média do CP do próprio sinal antes de qualquer análise. Uma forma simples de fazer isso e remover a tendência no sinal do CP é utilizar a função ‘detrend’ do Matlab 
[CPap=detrend(CPap); CPml=detrend(CPml)]. É possível aplicar um filtro passa-alta no sinal do CP. A escolha da frequência de corte desse 
filtro é bastante crítica e além do escopo deste texto. Uma vez que os dois procedimentos são realizados, diversas variáveis podem ser derivadas do sinal do CP.
//Deslocamento da oscilação total, DOT ‘Tamanho’ ou comprimento da trajetória do CP sobre a base de suporte;  DOT=sum(sqrt(CPap.^2+CPml.^2));
//Desvio-padrão, Dispersão do deslocamento do CP da posição média durante um intervalo de tempo; SDap=std(CPap); SDml=std(CPml);
//RMS (‘Root Mean Square’) mesmo resultado para RMS e desvio-padrão, se o sinal do CP tem média zero; RMSap=sqrt(sum(CPap.^2)/length(CPap)); RMSml=sqrt(sum(CPml.^2)/length(CPml);
//Amplitude de deslocamento do CP, Distância entre o deslocamento máximo e o mínimo do CP para cada direção; AdCPap=max(CPap) - min(CPap); AdCPml=max(CPml) - min(CPml);
//Velocidade média (VM), Determinação de quão rápidos foram os deslocamentos do CP VMap=sum(abs(diff(CPap)))*freq/length(CPap); VMml=sum(abs(diff(CPml)))*freq/length(CPml)
//Área [vec,val]=eig(cov(CPap,CPml)); Área=pi*prod(2.4478*sqrt(svd(val))) (a área é calculada a partir de método estatístico de análise dos componentes principais. Por meio dele, é possível o cálculo de uma elipse que engloba uma determinada porcentagem (por exemplo, 95%) dos dados do CP, sendo que os dois eixos da elipse são calculados a partir das medidas de dispersão dos sinais do CP).
//Velocidade média total (VMT) VMT=sum(sqrt(diff(CPap).^2+diff(CPml).^2))*freq/length(CPap) (A VMT é calculada a partir do deslocamento da oscilação total do CP nas duas direções dividido pelo tempo total da tentativa).
//Análise espectral: A análise de Fourier permite decompor um sinal qualquer como uma somatória de funções seno e coseno com diferentes amplitudes, frequências e fases. É possível obter informações sobre as frequências que compõem um sinal. Esse processo também é chamado de análise espectral, e o resultado dela é referido como o espectro do sinal original. Em termos práticos, a análise espectral é extremamente dependente do algoritmo e de seus parâmetros de entrada, o que dificulta a comparação dos resultados.
// Espectro: As frequências para um sinal do CP na direção ap e o código de Matlab são apresentadas a frente. A frequência predominante ou de pico é a que possui maior amplitude de todas as frequências que compõem o espectro. Baratto et al. sugerem que a banda de frequência com 80% da potência espectral é a que melhor caracteriza as alterações do sistema de controle postural. Além da análise nessas frequências, é comum a utilização da frequência média e frequência mediana do sinal. Para obter estimativas das características de frequência do sinal do CP, o método do periodograma de Welch pode ser utilizado no Matlab.
//Frequências: nfft=round(lenght(CPap)/2); [p,f]=psd(CPap,nfft,freq,nfft,round(nfft/2),'mean'); [m,peak]=max(p); area=cumtrapz(f,p); Find50=find(area>=.50*area(end)); Find80=find(area>=.80*area(end)); Fmedia=trapz(f,f,*p)/trapz(f,p); Fpico=f(peak); F50=f(Find50(1)); F80=f(Find80(1))
//Análise estrutural: decomposição do sinal de CP (movimento browniano, densidade de oscilação, não-estocástica)
//Collins e De Luca40 basearam-se na idéia de decomposição dos padrões de oscilação do sinal do CP em dois processos estocásticos modelados como ‘random walk’ ou movimento browniano: um processo de curta duração e um de longa duração. O movimento browniano é um processo estocástico em que, para cada instante de tempo, um passo é dado com amplitude fixa e direção randômica. Uma característica desse processo é que sua variância cresce linearmente com o tempo. Dessa forma, gráficos de difusão são construídos considerando-se
pares de amostras de dados do CP separados por um intervalo de tempo e computando a variância dos vetores correspondentes como uma função da amplitude do intervalo de tempo. Os autores propuseram que o controle postural humano seria composto por um circuito aberto
(que opera em intervalos de até cerca de 1 segundo) e um circuito fechado (que opera em intervalos maiores que cerca de 1 segundo).
//A análise estrutural proposta por Baratto et al. é baseada num conceito denominado curvas de densidade da oscilação. A idéia fundamental é que a estabilização postural é garantida pelo mecanismo de feedforward e, assim, o processo de controle é baseado em uma sequência de comandos motores antecipatórios. As curvas de densidade da oscilação são construídas pela contagem do número de amostras consecutivas da trajetória do CP que caem dentro de um círculo de raio conhecido. As curvas de densidade da oscilação são
caracterizadas por picos que representam instantes de tempo em que o momento de força no tornozelo e os comandos motores são relativamente estáveis e por vales que representam os instantes de tempo em que o momento de força no tornozelo muda rapidamente de um valor estável para outro. Duas variáveis seriam recomendadas na análise postural: a amplitude média dos picos e o tempo médio
entre os picos.
//A análise estrutural proposta por Duarte e Zatsiorsky é baseada na idéia de que a trajetória do CP não é puramente estocástica e que é possível identificar padrões consistentes por meio de uma análise no domínio espacial do estatocinesigrama e uma análise no domínio temporal do estabilograma. Tal análise é indicada para avaliação de tarefas de longa duração, em que é permitido ao indivíduo que está sendo avaliado realizar mudanças posturais, se assim desejar. Essas mudanças são geralmente observadas na postura natural, quando se está em pé e executando outra tarefa, por exemplo, conversando com alguém ou aguardando sua vez numa fila. Essa tarefa, quando investigada no laboratório, foi referida como postura irrestrita prolongada. Duarte e Zatsiorsky mostraram que quando o CP é apresen-
tado como uma série temporal, três padrões podem ser identificados: Shifting (tipo degrau): rápido deslocamento da posição média do CP de uma região à outra; Fidgeting (tipo pulso): rápido e grande deslocamento do CP e retorno para aproximadamente a mesma posição e Drifting (tipo rampa): contínuo e lento deslocamento da posição média do CP.
#calculate each variable for each trial
#variables, total área, total displacement, dimensions, elipse, dimensions on x and y (covariância), velocity in x and y . rms in x and r, load, pression and shear forces
#all trials, all test conditions
#new patient 
#all process
#indidual tab, line
#all group
#tab
#graphes, real data, on the base, statokinesiogram and stabilogram
#named compared groups
#estatística
#estatística, média, SP, grupo, p, gráficos, correlação, anova, big data, multiprofissional

//Reference: Marcos Duarte, Sandra M. S. F. Freitas. Revisão sobre posturografia baseada em plataforma de força para avaliação do equilíbrio. Rev Bras Fisioter, São Carlos, v. 14, n. 3, p. 183-92, maio/jun. 2010. ISSN 1413-3555
