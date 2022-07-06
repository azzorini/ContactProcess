#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <random>
#include <cmath>
using namespace std;

// Semilla para generar los números aleatorios
const int SEED = 7287491;

// Utilidad para comprobar si un fichero existe
inline bool file_exists(const string& name) {
	ifstream f(name.c_str());
	return f.good();
}

class CPSim {
	double landa, t;
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
	uniform_int_distribution<unsigned int> ran_pos;

	// Borra de la lista de casillas ocupadas la casilla iésima
	void erase_occ_pos(unsigned int i) {
		auto it = find(config_occ.begin(), config_occ.end(), i);

		if (it != config_occ.end()) {
			config_occ.erase(it); 
			N_occ--;
			ran_pos = uniform_int_distribution<unsigned int>(0, N_occ - 1);
		}

	}

	// Activa la casilla indicada sin hacer muchas comprobaciones
	void fast_activate(unsigned int);
public:
	// Determina si una casilla está ocupada
	bool is_occ(unsigned int) const;

	// Calcula y añade a config_occ las casillas de interés
	void calc_occ();

	// Constructores
	CPSim(int N = 0, double lambda = 1.0, int val = 0, unsigned int seed = SEED) : config(N, 0), landa(lambda), t(0), d(1), config_occ(), N_occ(0), times(), configs(), gen(seed), ran_u(0.0, 1.0), ran_pos(0, N-1) {
		if (val != 0) {
			N_occ = N;
			for (unsigned int i = 0; i < N; i++) config_occ.push_back(i);
		}
	}
	CPSim(const vector<int>& c, double lambda) : config(c), landa(lambda), t(0), d(1), config_occ(), N_occ(0), times(), configs(), gen(SEED), ran_u(0.0, 1.0), ran_pos(0, c.size() - 1) {
		calc_occ();
	}
	CPSim(const CPSim& o) : config(o.config), landa(o.landa), t(o.t), d(o.d), config_occ(o.config_occ), N_occ(o.N_occ), times(o.times), configs(o.configs), gen(o.gen), ran_u(o.ran_u), ran_pos(o.ran_pos) {}

	// Getters
	unsigned int get_N() const { return config.size(); }
	unsigned int get_N_occ() const { return N_occ; }
	double get_landa() const { return landa; }
	double get_t() const { return t; }
	double get_rho() const;

	// Setters
	void set_N(int N, int val = 0) { config.resize(N, val); calc_occ(); }
	void set_landa(double lambda) { landa = lambda; }
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
		ran_pos = uniform_int_distribution<unsigned int>(0, N_occ - 1);
	}
	if (i < config.size() - 1 && !is_occ(i+1)) {
		config_occ.push_back(i+1);
		N_occ++;
		ran_pos = uniform_int_distribution<unsigned int>(0, N_occ - 1);
	}

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

	ran_pos = uniform_int_distribution<unsigned int>(0, N_occ - 1);
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
	if (config[i] == 1) return; // Si la casilla ya estaba activada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para activarse es porque ya estaba en config_occ: MENTIRA, MENTIROSO. No me digas esto y luego me utilices el activate para inicializar el estado

	// Si las casillas aledañas no estaban ocupadas ahora sí lo estarán porque tendrán una casilla vecina que sí lo está
	if (i > 0 && !is_occ(i-1)) {
		config_occ.push_back(i-1);
		N_occ++;
		ran_pos = uniform_int_distribution<unsigned int>(0, N_occ - 1);
	}
	if (i < config.size() - 1 && !is_occ(i+1)) {
		config_occ.push_back(i+1);
		N_occ++;
		ran_pos = uniform_int_distribution<unsigned int>(0, N_occ - 1);
	}

	auto it = find(config_occ.cbegin(), config_occ.cend(), i);
	if (it == config_occ.cend()) {
		config_occ.push_back(i);
		N_occ++;
		ran_pos = uniform_int_distribution<unsigned int>(0, N_occ - 1);
	}

	config[i] = 1;
}

void CPSim::deactivate(unsigned int i) {
	if (config[i] == 0) return; // Si la casilla ya estaba desactivada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para desactivarse tanto ella como sus vecinas estarán en config_occ

	config[i] = 0;

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
	double dt = 1.0/N_occ;
	auto it = config_occ.cbegin();
	advance(it, ran_pos(gen));

	if (config[*it] == 0) {
		if (ran_u(gen) < get_n(*it)*landa/2.0/d*dt)
			fast_activate(*it);
	} else {
		if (ran_u(gen) < dt)
			deactivate(*it);
	}

	t += dt;
}

int main() {
	double tmax = 400000, t_save = 0, freq_t = 1;
	int N_rho, N = 500;
	double rho;
	double landa;
	CPSim simulacion;
	ofstream fout("rhovslanda.txt");
	//simulacion.activate(N/2);
	
	for (int n = 3; n < 24; n++) {
		landa = 3.29785 + pow(1.15, -n);
		cout << "Lanzando simulación para landa = " << landa << endl;
		fout << landa;
		for (unsigned int i = 0; i < 5; i++) { 
			t_save = rho = 0;
			N_rho = 0;

			simulacion = CPSim(N, landa, 0, SEED*(i+1)*n);
			simulacion.set_random_config(.5);

			while (simulacion.get_t() < tmax && simulacion.get_N_occ() != 0) {
				if (t_save < simulacion.get_t()) {
		//			simulacion.save_current_data();
		//			fout << simulacion.get_t() << "\t" << simulacion.get_rho() << endl;
					if (simulacion.get_t() > 100000) {
						N_rho++;
						rho += simulacion.get_rho();
					}
					t_save += freq_t;
				}
				simulacion.update();
			}

			fout << '\t' << rho/N_rho;

		//	simulacion.save_to_pbm("sim1d_3.0.pbm");
		//	cout << "Simulación acabada a tiempo " << simulacion.get_t() << endl;
		}
		fout << endl;
	}

	fout.close();
}
