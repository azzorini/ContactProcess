#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <random>
#include <cmath>
#include <cfloat>
#include <thread>
#include <syncstream>
using namespace std;

// Semilla para generar los números aleatorios
const int SEED = 672239;

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
	unsigned int get_data_N() const { return times.size(); }
	double get_data_t(unsigned int i) const { return times[i]; }
	int get_data(unsigned int i, unsigned int j) const { return configs[i][j]; }

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
	if (i > 0) {
		if (!is_occ(i-1)) {
			config_occ.push_back(i-1);
			N_occ++;
		}
	} else if (!is_occ(config.size()-1)) {
		config_occ.push_back(config.size()-1);
		N_occ++;
	}
	if (i < config.size() - 1) {
		if (!is_occ(i+1)) {
			config_occ.push_back(i+1);
			N_occ++;
		}
	} else if (!is_occ(0)) {
		config_occ.push_back(0);
		N_occ++;
	}

	W += 1 - get_n(i)*landa/2/d;
	if (i > 0) {
		if (config[i-1] == 0) W += landa/2/d;
	} else {
		if (config[config.size()-1] == 0) W += landa/2/d;
	}

	if (i < config.size() - 1) {
		if (config[i+1] == 0) W += landa/2/d;
	} else {
		if (config[0] == 0) W += landa/2/d;
	}

	config[i] = 1;
}

bool CPSim::is_occ(unsigned int i) const {
	if (config[i] != 0) return true;
	else {
		if (i > 0) {
			if (config[i-1] != 0) return true;
		} else {
			if (config.back() != 0) return true;
		}

		if (i < config.size() - 1) {
			if (config[i+1] != 0) return true;
		} else {
			if (config.front() != 0) return true;
		}
	}

	return false;
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
	for (auto it = config_occ.cbegin(); it != config_occ.cend(); it++) r += config[*it];
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
	else n += config.back();

	if (i < config.size() - 1) n += config[i+1];
	else n += config.front();

	return n;
}

void CPSim::activate(unsigned int i) {
	bool added = false;

	if (config[i] == 1) return; // Si la casilla ya estaba activada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para activarse es porque ya estaba en config_occ: MENTIRA, MENTIROSO. No me digas esto y luego me utilices el activate para inicializar el estado

	// Si las casillas aledañas no estaban ocupadas ahora sí lo estarán porque tendrán una casilla vecina que sí lo está
	if (i > 0) {
		if (!is_occ(i-1)) {
			config_occ.push_back(i-1);
			N_occ++;
		}
	} else if (!is_occ(config.size()-1)) {
		config_occ.push_back(config.size()-1);
		N_occ++;
	}
	if (i < config.size() - 1) {
		if (!is_occ(i+1)) {
			config_occ.push_back(i+1);
			N_occ++;
		}
	} else if (!is_occ(0)) {
		config_occ.push_back(0);
		N_occ++;
	}

	auto it = find(config_occ.cbegin(), config_occ.cend(), i);
	if (it == config_occ.cend()) {
		config_occ.push_back(i);
		added = true;
		N_occ++;
	}

	W += (added ? 1 : 1 - get_n(i)*landa/2/d);
	if (i > 0) {
		if (config[i-1] == 0) W += landa/2/d;
	} else {
		if (config[config.size()-1] == 0) W += landa/2/d;
	}

	if (i < config.size() - 1) {
		if (config[i+1] == 0) W += landa/2/d;
	} else {
		if (config[0] == 0) W += landa/2/d;
	}

	config[i] = 1;
}

void CPSim::deactivate(unsigned int i) {
	if (config[i] == 0) return; // Si la casilla ya estaba desactivada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para desactivarse tanto ella como sus vecinas estarán en config_occ

	config[i] = 0;

	W += get_n(i)*landa/2/d - 1;
	if (i > 0) {
		if (config[i-1] == 0) W -= landa/2/d;
	} else {
		if (config.back() == 0) W -= landa/2/d;
	}

	if (i < config.size() - 1) {
		if (config[i+1] == 0) W -= landa/2/d;
	} else {
		if (config.front() == 0) W -= landa/2/d;
	}

	if (!is_occ(i)) {
		erase_occ_pos(i);
	}
	if (i > 0) {
		if (!is_occ(i-1)) erase_occ_pos(i-1);
	} else {
		if (!is_occ(config.size()-1)) erase_occ_pos(config.size()-1);
	}

	if (i < config.size() - 1) {
		if (!is_occ(i+1)) erase_occ_pos(i+1);
	} else {
		if (!is_occ(0)) erase_occ_pos(0);
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

// Función que nos devuelve la correlación espacial de una isla en una simulación por encima del
// landa crítico
unsigned int get_spatial_correlation_length(const CPSim& sim, const vector< vector<int> >& data, unsigned int n_island) {
	unsigned int N = sim.get_N(), N_data = sim.get_data_N();
	unsigned int i, j;
	int hole_length, max_length = 0;
	bool in_hole = false;

	// Bucle por todos los pasos temporales
	for (i = 0; i < N_data; i++) {
		for (j = 0; j < N; j++) {
			if (in_hole) {
				if (data[i][j] == n_island) hole_length++;
				else {
					in_hole = false;
					if (hole_length > max_length) max_length = hole_length;
				}
			} else {
				if (data[i][j] == n_island) {
					in_hole = true;
					hole_length = 1;
				}
			}
		}

		if (in_hole) {
			j = 0;
			while (data[i][j] == n_island) {
				hole_length++; j++;
			}
			in_hole = false;
			if (hole_length > max_length) max_length = hole_length;
		}
	}

	return max_length;
}


// Función que nos devuelve la correlación temporal de una isla en una simulación por encima del
// landa crítico
double get_time_correlation_length(const CPSim& sim, const vector< vector<int> >& data, unsigned int n_island) {
	unsigned int N = sim.get_N(), N_data = sim.get_data_N();
	double first_t, max_length = 0;
	bool in_hole = false;

	// Bucle por todos los pasos temporales
	for (unsigned int j = 0; j < N; j++) {
		for (unsigned int i = 0; i < N_data; i++) {
			if (in_hole) {
				if (data[i][j] != n_island) {
					in_hole = false;
					if (sim.get_data_t(i) - first_t > max_length) max_length = sim.get_data_t(i) - first_t;
				}
			} else {
				if (data[i][j] == n_island) {
					in_hole = true;
					first_t = sim.get_data_t(i);
				}
			}
		}

		if (in_hole) {
			in_hole = false;
			if (sim.get_data_t(N_data - 1) - first_t > max_length) max_length = sim.get_data_t(N_data - 1) - first_t;
		}
	}

	return max_length;
}

struct Position {
	unsigned int i, j;

	// Constructor
	Position(unsigned int x = 0, unsigned int y = 0) : i(x), j(y) {}

	// Operador comparación: Lo necesitaremos para poder usar std::find
	bool operator==(const Position& p) const { return i == p.i && j == p.j; }
};

void add_island(vector< vector<int> >& data, double& size, const CPSim& sim, unsigned int n_island, unsigned int i, unsigned int j) {
	vector<Position> isla;
	Position aux;

	if (data[i][j] != -1) return; // Esta casilla ya ha sido contabilizada

	if (sim.get_data(i, j) == 0) {
		data[i][j] = n_island;
		isla.push_back(Position(i, j));
		size += sim.get_data_t(i+1) - sim.get_data_t(i);
	} else data[i][j] = 0;

	while (!isla.empty()) {
		aux = isla.back();
		isla.pop_back();

		if (aux.i > 0 && data[aux.i-1][aux.j] == -1) {
			if (sim.get_data(aux.i-1, aux.j) == 0) {
				data[aux.i-1][aux.j] = n_island;
				isla.push_back(Position(aux.i-1, aux.j));
				size += sim.get_data_t(aux.i) - sim.get_data_t(aux.i-1);

			} else data[aux.i-1][aux.j] = 0;
		} 

		if (aux.i < sim.get_data_N() - 2 && data[aux.i+1][aux.j] == -1) {
			if (sim.get_data(aux.i+1, aux.j) == 0) {
				data[aux.i+1][aux.j] = n_island;
				isla.push_back(Position(aux.i+1, aux.j));
				size += sim.get_data_t(aux.i+2) - sim.get_data_t(aux.i+1);
			} else data[aux.i+1][aux.j] = 0;
		}

		if (aux.j > 0) {
			if (data[aux.i][aux.j-1] == -1) {
				if (sim.get_data(aux.i, aux.j-1) == 0) {
					data[aux.i][aux.j-1] = n_island;
					isla.push_back(Position(aux.i, aux.j-1));
					size += sim.get_data_t(aux.i+1) - sim.get_data_t(aux.i);
				} else data[aux.i][aux.j-1] = 0;
			}
		} else {
			if (data[aux.i][sim.get_N() - 1] == -1) {
				if (sim.get_data(aux.i, sim.get_N() - 1) == 0) {
					data[aux.i][sim.get_N() - 1] = n_island;
					isla.push_back(Position(aux.i, sim.get_N() - 1));
					size += sim.get_data_t(aux.i+1) - sim.get_data_t(aux.i);
				} else data[aux.i][sim.get_N() - 1] = 0;
			}
		}

		if (aux.j < sim.get_N() - 1) {
			if (data[aux.i][aux.j+1] == -1) {
				if (sim.get_data(aux.i, aux.j+1) == 0) {
					data[aux.i][aux.j+1] = n_island;
					isla.push_back(Position(aux.i, aux.j+1));
					size += sim.get_data_t(aux.i+1) - sim.get_data_t(aux.i);
				} else data[aux.i][aux.j+1] = 0;
			}
		} else {
			if (data[aux.i][0] == -1) {
				if (sim.get_data(aux.i, 0) == 0) {
					data[aux.i][0] = n_island;
					isla.push_back(Position(aux.i, 0));
					size += sim.get_data_t(aux.i+1) - sim.get_data_t(aux.i);
				} else data[aux.i][0] = 0;
			}
		}
	}

	// Versión recursiva
	//if (sim.get_data(i, j) == 0) {
		//data[i][j] = n_island;
		//size += sim.get_data_t(i+1) - sim.get_data_t(i);
//
		//if (i > 0) add_island(data, size, sim, n_island, i-1, j);
//
		//if (i < sim.get_data_N() - 2) add_island(data, size, sim, n_island, i+1, j);
//
		//if (j > 0) add_island(data, size, sim, n_island, i, j-1);
		//else add_island(data, size, sim, n_island, i, sim.get_N() - 2);
//
		//if (j < sim.get_N() - 1) add_island(data, size, sim, n_island, i, j+1);
		//else add_island(data, size, sim, n_island, i, 0);
	//} else data[i][j] = 0;

	return;
}

// Función que calcula la isla de mayor tamaño y sus dimensiones
void get_corrleation_lengths(unsigned int& space_cl, double& time_cl, const CPSim& sim) {
	double size, max_size = 0;
	unsigned int island_counter = 0, max_island = 0;
	vector< vector<int> > data(sim.get_data_N(), vector<int>(sim.get_N(), -1));

	// Ponemos en data todas las islas y nos guardamos cual es la isla de mayor tamaño
	for (unsigned int i = 0; i < sim.get_data_N() - 1; i++) {
		for (unsigned int j = 0; j < sim.get_N(); j++) {
			if (data[i][j] == -1) {
				if (sim.get_data(i, j) != 0) data[i][j] = 0;
				else {
					island_counter++; size = .0;
					add_island(data, size, sim, island_counter, i, j);
					if (size > max_size) {
						max_size = size;
						max_island = island_counter;
					}
				}
			}
		}
	}

	space_cl = get_spatial_correlation_length(sim, data, max_island);
	time_cl = get_time_correlation_length(sim, data, max_island);
}

void simula_landas(int N, double tmax, int start, int stop, const string& filename) {
	double t_save, freq_t = .5;
	double landa, rho_stat;
	unsigned int space_cl; double time_cl;
	CPSim simulacion;
	ofstream ftime(("PBC_TCL_m3_1I_" + filename).c_str()), fspace(("PBC_SCL_m3_1I_" + filename).c_str());
	
	for (int n = start; n < stop; n++) {
		landa = 3.29785 + pow(1.15, -n);
		rho_stat = exp(-0.406)*pow(pow(1.15, -n), 0.274); // Estimación sacada del ajuste para la densidad
		{
			osyncstream out{cout};
			out << "Lanzando simulaciones para landa = " << landa << endl;
		}
		ftime << landa; fspace << landa;
		for (unsigned int i = 0; i < 5; i++) { 
			t_save = 0;

			simulacion = CPSim(N, landa, 0, SEED*(i+1)*(n+20));
			simulacion.set_random_config(rho_stat);

			while (simulacion.get_t() < tmax && simulacion.get_N_occ() != 0) {
				if (t_save < simulacion.get_t()) {
					simulacion.save_current_data();
					t_save += freq_t;
				}
				simulacion.update();
			}

			get_corrleation_lengths(space_cl, time_cl, simulacion);
			fspace << '\t' << space_cl;
			ftime << '\t' << time_cl;
		}
		ftime << endl; fspace << endl;
	}

	ftime.close(); fspace.close();

}

int main() {
	int N = 1000;
	double tmax = 40000;

	jthread t1(simula_landas, N, tmax, -16, -3, "landas_grandes.txt");
	jthread t2(simula_landas, N, tmax, -3, 11, "landas_medianos.txt");
	jthread t3(simula_landas, N, tmax, 11, 25, "landas_pequenyos.txt");
}
