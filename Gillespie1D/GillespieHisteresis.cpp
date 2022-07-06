#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <random>
#include <cmath>
#include <cfloat>
using namespace std;

// Semilla para generar los números aleatorios
const int SEED = 123456;

// Utilidad para comprobar si un fichero existe
inline bool file_exists(const string& name) {
	ifstream f(name.c_str());
	return f.good();
}

class CPSim {
	double landa, t, W;
	int d;
	vector<int> config;

	// Usaremos estos objetos para ir almacenando los datos que queramos sacar
	vector<double> times;
	vector< vector<int> > configs;

	list<unsigned int> config_occ;
	unsigned int N_occ;

	// Usaremos esto para generar los números aleatorios
	mt19937 gen;
	uniform_real_distribution<double> ran_u;

	// Borra de la lista de casillas ocupadas la casilla iésima
	void erase_occ_pos(unsigned int i) {
		auto it = find(config_occ.begin(), config_occ.end(), i);

		if (it != config_occ.end()) {
			config_occ.erase(it); 
			N_occ--;
		}

	}

	// Activa la casilla indicada sin hacer muchas comprobaciones
	void fast_activate(unsigned int);
public:
	// Determina si una casilla está ocupada
	bool is_occ(unsigned int) const;

	// Calcula y añade a config_occ las casillas de interés
	void calc_occ();

	// Calcula la probabilidad de cambiar de estado
	double calc_W() const;

	// Constructores
	CPSim(int N = 0, double lambda = 1.0, int val = 0, unsigned int seed = SEED) : config(N, val), landa(lambda), t(0), W(0), d(1), config_occ(), N_occ(0), times(), configs(), gen(seed), ran_u(DBL_TRUE_MIN, 1.0) {
		if (val != 0) {
			N_occ = W = N;
			for (unsigned int i = 0; i < N; i++) config_occ.push_back(i);
		}
	}
	CPSim(const vector<int>& c, double lambda) : config(c), landa(lambda), t(0), W(0), d(1), config_occ(), N_occ(0), times(), configs(), gen(SEED), ran_u(DBL_TRUE_MIN, 1.0) {
		calc_occ();
		W = calc_W();
	}
	CPSim(const CPSim& o) : config(o.config), landa(o.landa), t(o.t), W(o.W), d(o.d), config_occ(o.config_occ), N_occ(o.N_occ), times(o.times), configs(o.configs), gen(o.gen), ran_u(o.ran_u) {}

	// Getters
	unsigned int get_N() const { return config.size(); }
	unsigned int get_N_occ() const { return N_occ; }
	double get_landa() const { return landa; }
	double get_t() const { return t; }
	double get_W() const { return W; }
	double get_rho() const;

	// Setters
	void set_N(int N, int val = 0) { config.resize(N, val); calc_occ(); }
	void set_landa(double lambda) { landa = lambda; W = calc_W(); }
	void set_t(double time) { t = time; }
	// Inicializa a una configuración aleatoria siguendo aproximadamente una cierta densidad
	// si esta densidad es mayor que 1 esta configuración será todos activos
	// en caso de que sea menor que 0 la configuración será todos inactivos
	void set_random_config(double);

	// Operador []
	const int& operator[](unsigned int i) const { return config[i]; }

	// Guardamos los datos en el tiempo actual para sacarlos a posteriori
	void save_current_data() {
		times.push_back(t);
		configs.push_back(config);
	}

	// Limpiamos todos los datos que teníamos almacenandos para sacar
	void clear_data() {
		times.clear();
		config.clear();
	}

	// Guardamos el fichero en formato pbm
	void save_to_pbm(const string&) const;

	// Guardamos el fichero en formato txt
	void save_to_txt(const string&) const;

	// Obtiene el número de vecinos infectados del sitio i
	int get_n(unsigned int) const;

	// Activa la casilla indicada y actualiza la información que sea necesaria
	void activate(unsigned int);

	// Descativa la casilla indicada y actualiza la información que sea necesaria
	void deactivate(unsigned int);

	// Actualizamos el sistema de t a t + dt
	void update();
};

void CPSim::fast_activate(unsigned int i) {
	if (config[i] == 1) return; // Si la casilla ya estaba activada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para activarse es porque ya estaba en config_occ: MENTIRA, MENTIROSO. Pero en esta función privada así se asumirá

	// Si las casillas aledañas no estaban ocupadas ahora sí lo estarán porque tendrán una casilla vecina que sí lo está
	if (i > 0 && !is_occ(i-1)) {
		config_occ.push_back(i-1);
		N_occ++;
	}
	if (i < config.size() - 1 && !is_occ(i+1)) {
		config_occ.push_back(i+1);
		N_occ++;
	}

	W += 1 - get_n(i)*landa/2/d;
	if (i > 0 && config[i-1] == 0) W += landa/2/d;
	if (i < config.size() - 1 && config[i+1] == 0) W += landa/2/d;

	config[i] = 1;
}

bool CPSim::is_occ(unsigned int i) const {
	if (config[i] != 0) {
		return true;
	} else if (i > 0 && config[i-1] != 0) {
		return true;
	} else if (i < config.size() - 1 && config[i+1] != 0) {
		return true;
	} else return false;
}

void CPSim::calc_occ() {
	config_occ.clear();
	N_occ = 0;

	for (unsigned int i = 0; i < config.size(); i++) {
		if (is_occ(i)) {
			config_occ.push_back(i);
			N_occ++;
		}
	}

}

double CPSim::calc_W() const {
	double r = .0;

	for (auto it = config_occ.cbegin(); it != config_occ.cend(); it++) {
		if (config[*it] == 0) r += get_n(*it)*landa/2/d;
		else r += 1;
	}

	return r;
}

double CPSim::get_rho() const {
	unsigned int r = 0;
	for (unsigned int i = 0; i < config.size(); i++) r += config[i];
	return ((double) r)/config.size();
}

void CPSim::set_random_config(double rho = .5) {
	for (unsigned int i = 0; i < config.size(); i++) {
		if (ran_u(gen) < rho) activate(i);
		else deactivate(i);
	}
}

void CPSim::save_to_pbm(const string& filename = "output.pbm") const {
	string line, contenido_previo = "";
	int prev_times = 0;
	ofstream fout; ifstream fin;

	if (file_exists(filename)) {
		fin.open(filename.c_str());
		
		getline(fin, line); getline(fin, line); // Pasamos de las dos primeras líneas (P1 + comentario)
		fin >> prev_times >> prev_times; getline(fin, line); // Leemos el número de tiempos del fichero anterior y obviamos el resto

		// Leemos todos los datos del fichero anterior
		while (getline(fin, line))
			contenido_previo += line + "\n";

		fin.close();
	}

	fout.open(filename.c_str());

	fout << "P1\n# Simulación con landa = " << landa
	     << " y N = " << config.size() << endl
	     << config.size() << " " << configs.size() + prev_times << endl
	     << contenido_previo;

	for (unsigned int i = 0; i < configs.size(); i++) {
		fout << configs[i][0];
		for (unsigned int j = 1; j < configs[i].size(); j++)
			fout << " " << configs[i][j];
		fout << endl;
	}

	fout.close();
}

void CPSim::save_to_txt(const string& filename = "output.txt") const {
	bool file_created = file_exists(filename);
	fstream fout(filename.c_str(), fstream::app | fstream::out);

	// Si el archivo no existía ya añadimos el comentario al principio
	if (!file_created)
		fout << "# t\tx_0\tx_1\t...\tx_n\n";

	for (unsigned int i = 0; i < configs.size(); i++) {
		fout << configs[i][0];
		for (unsigned int j = 1; j < configs[i].size(); j++)
			fout << "\t" << configs[i][j];
		fout << endl;
	}

	fout.close();
}

int CPSim::get_n(unsigned int i) const {
	int n = 0;

	// Sumamos los vecinos infectados teniendo cuidado con los bordes
	if (i > 0) n += config[i-1];
	if (i < config.size() - 1) n += config[i+1];

	return n;
}

void CPSim::activate(unsigned int i) {
	bool added = false;

	if (config[i] == 1) return; // Si la casilla ya estaba activada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para activarse es porque ya estaba en config_occ: MENTIRA, MENTIROSO. No me digas esto y luego me utilices el activate para inicializar el estado

	// Si las casillas aledañas no estaban ocupadas ahora sí lo estarán porque tendrán una casilla vecina que sí lo está
	if (i > 0 && !is_occ(i-1)) {
		config_occ.push_back(i-1);
		N_occ++;
	}
	if (i < config.size() - 1 && !is_occ(i+1)) {
		config_occ.push_back(i+1);
		N_occ++;
	}

	auto it = find(config_occ.cbegin(), config_occ.cend(), i);
	if (it == config_occ.cend()) {
		config_occ.push_back(i);
		added = true;
		N_occ++;
	}

	W += (added ? 1 : 1 - get_n(i)*landa/2/d);
	if (i > 0 && config[i-1] == 0) W += landa/2/d;
	if (i < config.size() - 1 && config[i+1] == 0) W += landa/2/d;

	config[i] = 1;
}

void CPSim::deactivate(unsigned int i) {
	if (config[i] == 0) return; // Si la casilla ya estaba desactivada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para desactivarse tanto ella como sus vecinas estarán en config_occ

	config[i] = 0;

	W += get_n(i)*landa/2/d - 1;
	if (i > 0 && config[i-1] == 0) W -= landa/2/d;
	if (i < config.size() - 1 && config[i+1] == 0) W -= landa/2/d;

	if (!is_occ(i)) {
		erase_occ_pos(i);
	}
	if (i > 0 && !is_occ(i-1)) {
		erase_occ_pos(i-1);
	}
	if (i < config.size() - 1 && !is_occ(i+1)) {
		erase_occ_pos(i+1);
	}
}

void CPSim::update() {
	if (N_occ == 0) return;

	double p, r;
	t += -log(ran_u(gen))/W;

	auto it = config_occ.cbegin();
	p = ran_u(gen)*W;
	r = (config[*it] == 0 ? get_n(*it)*landa/2/d : 1);
	while (r < p && ++it != config_occ.cend()) {
		r += (config[*it] == 0 ? get_n(*it)*landa/2/d : 1);
	}

	if (it != config_occ.cend()) {
		if (config[*it] == 0) fast_activate(*it);
		else deactivate(*it);
	}
}

int main() {
	int N = 1000;
	unsigned int n_it = 6000;
	unsigned int n_it_per_lambda = 50;
	double landa = 3.29785 - 1; // Vamos a ir desde 1 por debajo del punto crítico hasta 5 por arriba
	double d_lambda = .001;
	double rho;

	CPSim simulacion;

	simulacion = CPSim(N, landa, 0, SEED);
	simulacion.activate(N/2);

	ofstream fout("ciclo_histeresis_50.txt");
	fout << "#Lambda\tRho\n";

	// Amos parriba
	for (unsigned int i = 0; i < n_it; i++) {
		rho = .0;
		for (unsigned int j = 0; j < n_it_per_lambda; j++) {
			simulacion.update();
			if (simulacion.get_N_occ() == 0) simulacion.activate(N/2);
			rho += simulacion.get_rho();
		}
		fout << simulacion.get_landa() << '\t' << rho/n_it_per_lambda << endl;
		landa += d_lambda;
		simulacion.set_landa(landa);
	}

	// Amos pabajo
	for (unsigned int i = 0; i < n_it; i++) {
		rho = .0;
		for (unsigned int j = 0; j < n_it_per_lambda; j++) {
			simulacion.update();
			if (simulacion.get_N_occ() == 0) simulacion.activate(N/2);
			rho += simulacion.get_rho();
		}
		fout << simulacion.get_landa() << '\t' << rho/n_it_per_lambda << endl;
		landa -= d_lambda;
		simulacion.set_landa(landa);
	}

	fout.close();
}
