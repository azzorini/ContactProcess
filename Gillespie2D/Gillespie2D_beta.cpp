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
const int SEED = 546;

// Utilidad para comprobar si un fichero existe
inline bool file_exists(const string& name) {
	ifstream f(name.c_str());
	return f.good();
}

struct Position {
	unsigned int i, j;

	// Constructor
	Position(unsigned int x = 0, unsigned int y = 0) : i(x), j(y) {}

	// Operador comparación: Lo necesitaremos para poder usar std::find
	bool operator==(const Position& p) const { return i == p.i && j == p.j; }
};

class CPSim {
	double landa, t, W;
	int d;
	vector< vector<int> > config;

	// Usaremos estos objetos para ir almacenando los datos que queramos sacar
	vector<double> times;
	vector< vector< vector<int> > > configs;

	list<Position> config_occ;
	unsigned int N_occ;

	// Usaremos esto para generar los números aleatorios
	mt19937 gen;
	uniform_real_distribution<double> ran_u;

	// Borra de la lista de casillas ocupadas la casilla iésima
	void erase_occ_pos(const Position& p) {
		auto it = find(config_occ.begin(), config_occ.end(), p);

		if (it != config_occ.end()) {
			config_occ.erase(it); 
			N_occ--;
		}

	}

	// Activa la casilla indicada sin hacer muchas comprobaciones
	void fast_activate(const Position&);
public:
	// Determina si una casilla está ocupada
	bool is_occ(const Position&) const;

	// Calcula y añade a config_occ las casillas de interés
	void calc_occ();

	// Calcula la probabilidad de cambiar de estado
	double calc_W() const;

	// Constructores
	CPSim(int N = 0, double lambda = 1.0, int val = 0, unsigned int seed = SEED) : config(N, vector<int>(N, val)), landa(lambda), t(0), W(0), d(2), config_occ(), N_occ(0), times(), configs(), gen(seed), ran_u(DBL_TRUE_MIN, 1.0) {
		if (val != 0) {
			N_occ = W = N*N;
			for (unsigned int i = 0; i < N; i++)
				for (unsigned int j = 0; j < N; j++) config_occ.push_back(Position(i, j));
		}
	}
	CPSim(const vector< vector<int> >& c, double lambda, unsigned int seed = SEED) : config(c), landa(lambda), t(0), W(0), d(1), config_occ(), N_occ(0), times(), configs(), gen(seed), ran_u(DBL_TRUE_MIN, 1.0) {
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
	void set_N(int, int);
	void set_landa(double lambda) { landa = lambda; W = calc_W(); }
	void set_t(double time) { t = time; }
	// Inicializa a una configuración aleatoria siguendo aproximadamente una cierta densidad
	// si esta densidad es mayor que 1 esta configuración será todos activos
	// en caso de que sea menor que 0 la configuración será todos inactivos
	void set_random_config(double);

	// Operador []
	const int& operator[](const Position& p) const { return config[p.i][p.j]; }

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

	// Guardamos la proyección en uno de los ejes en un fichero en formato pbm
	void save_to_pbm(const string&, unsigned int) const;

	// Guardamos el fichero en formato txt
	void save_to_txt(const string&) const;

	// Guardamos la proyección en uno de los ejes en un fichero en formato txt
	void save_to_txt_proyected(const string&, unsigned int) const;

	// Obtiene el número de vecinos infectados del sitio i
	int get_n(const Position&) const;

	// Activa la casilla indicada y actualiza la información que sea necesaria
	void activate(const Position&);

	// Descativa la casilla indicada y actualiza la información que sea necesaria
	void deactivate(const Position&);

	// Actualizamos el sistema de t a t + dt
	void update();
};

void CPSim::fast_activate(const Position& p) {
	if (config[p.i][p.j] == 1) return; // Si la casilla ya estaba activada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para activarse es porque ya estaba en config_occ: MENTIRA, MENTIROSO. Pero en esta función privada así se asumirá

	// Si las casillas aledañas no estaban ocupadas ahora sí lo estarán porque tendrán una casilla vecina que sí lo está
	if (p.i > 0 && !is_occ(Position(p.i - 1, p.j))) {
		config_occ.push_back(Position(p.i - 1, p.j));
		N_occ++;
	}
	if (p.i < config.size() - 1 && !is_occ(Position(p.i + 1, p.j))) {
		config_occ.push_back(Position(p.i + 1, p.j));
		N_occ++;
	}
	if (p.j > 0 && !is_occ(Position(p.i, p.j - 1))) {
		config_occ.push_back(Position(p.i, p.j - 1));
		N_occ++;
	}
	if (p.j < config[0].size() - 1 && !is_occ(Position(p.i, p.j + 1))) {
		config_occ.push_back(Position(p.i, p.j + 1));
		N_occ++;
	}

	W += 1 - get_n(p)*landa/2/d;
	if (p.i > 0 && config[p.i-1][p.j] == 0) W += landa/2/d;
	if (p.i < config.size() - 1 && config[p.i+1][p.j] == 0) W += landa/2/d;
	if (p.j > 0 && config[p.i][p.j-1] == 0) W += landa/2/d;
	if (p.j < config.size() - 1 && config[p.i][p.j+1] == 0) W += landa/2/d;

	config[p.i][p.j] = 1;
}

bool CPSim::is_occ(const Position& p) const {
	if (config[p.i][p.j] != 0) {
		return true;
	} else if (p.i > 0 && config[p.i-1][p.j] != 0) {
		return true;
	} else if (p.i < config.size() - 1 && config[p.i+1][p.j] != 0) {
		return true;
	} else if (p.j > 0 && config[p.i][p.j-1] != 0) {
		return true;
	} else if (p.j < config.size() - 1 && config[p.i][p.j+1] != 0) {
		return true;
	} else return false;
}

void CPSim::calc_occ() {
	Position p;
	config_occ.clear();
	N_occ = 0;

	for (unsigned int i = 0; i < config.size(); i++) {
		p.i = i;
		for (unsigned int j = 0; j < config[i].size(); j++) {
			p.j = j;
			if (is_occ(p)) {
				config_occ.push_back(p);
				N_occ++;
			}
		}
	}

}

double CPSim::calc_W() const {
	double r = .0;

	for (auto it = config_occ.cbegin(); it != config_occ.cend(); it++) {
		if (config[it->i][it->j] == 0) r += get_n(*it)*landa/2/d;
		else r += 1;
	}

	return r;
}

double CPSim::get_rho() const {
	unsigned int r = 0;
	for (unsigned int i = 0; i < config.size(); i++)
		for (unsigned int j = 0; j < config[i].size(); j++) r += config[i][j];
	return ((double) r)/config.size()/config.size();
}


void CPSim::set_N(int N, int val = 0) {
	int it_limit, ult_N = config.size();
	config.resize(N, vector<int>(N, val));
	calc_occ();
	W = calc_W();
	it_limit = (N > ult_N ? ult_N : N);
	for (unsigned int i = 0; i < it_limit; i++)
		config[i].resize(N, val);
}

void CPSim::set_random_config(double rho = .5) {
	for (unsigned int i = 0; i < config.size(); i++)
		for (unsigned int j = 0; j < config[i].size(); j++) {
			if (ran_u(gen) < rho) activate(Position(i, j));
			else deactivate(Position(i, j));
	}
}

void CPSim::save_to_pbm(const string& filename = "output.pbm", unsigned int ax = 0) const {
	string line, contenido_previo = "";
	int prev_times = 0, N = config.size();
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
	     << " y N = " << N << endl
	     << N << " " << configs.size() + prev_times << endl
	     << contenido_previo;

	for (unsigned int i = 0; i < configs.size(); i++) {
		if (ax == 0) fout << configs[i][0][N/2];
		else fout << configs[i][N/2][0];

		for (unsigned int j = 1; j < N; j++) {
			if (ax == 0) fout << ' ' << configs[i][j][N/2];
			else fout << ' ' << configs[i][N/2][j];
		}
		fout << endl;
	}

	fout.close();
}

void CPSim::save_to_txt(const string& filename = "output.txt") const {
	bool file_created = file_exists(filename);
	fstream fout(filename.c_str(), fstream::app | fstream::out);

	// Si el archivo no existía ya añadimos el comentario al principio
	if (!file_created)
		fout << "# x\ty\tt\tval\n";

	for (unsigned int i = 0; i < configs.size(); i++)
		for (unsigned int j = 0; j < configs[i].size(); j++)
			for (unsigned int k = 0; k < configs[i][j].size(); k++)
				fout << j << '\t'  << k << '\t' << times[i]<< '\t' << configs[i][j][k] << '\n';
	

	fout.close();
}

void CPSim::save_to_txt_proyected(const string& filename = "output.txt", unsigned int ax = 0) const {
	unsigned int N = config.size();
	bool file_created = file_exists(filename);
	fstream fout(filename.c_str(), fstream::app | fstream::out);

	// Si el archivo no existía ya añadimos el comentario al principio
	if (!file_created)
		fout << "# t\tx_0\tx_1\t...\tx_n\n";

	for (unsigned int i = 0; i < configs.size(); i++) {
		if (ax == 0) fout << configs[i][0][N/2];
		else fout << configs[i][N/2][0];

		for (unsigned int j = 1; j < N; j++) {
			if (ax == 0) fout << '\t' << configs[i][j][N/2];
			else fout << '\t' << configs[i][N/2][j];
		}
		fout << endl;
	}

	fout.close();
}

int CPSim::get_n(const Position& p) const {
	int n = 0;

	// Sumamos los vecinos infectados teniendo cuidado con los bordes
	if (p.i > 0) n += config[p.i-1][p.j];
	if (p.i < config.size() - 1) n += config[p.i+1][p.j];
	if (p.j > 0) n += config[p.i][p.j-1];
	if (p.j < config.size() - 1) n += config[p.i][p.j+1];

	return n;
}

void CPSim::activate(const Position& p) {
	bool added = false;

	if (config[p.i][p.j] == 1) return; // Si la casilla ya estaba activada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para activarse es porque ya estaba en config_occ: MENTIRA, MENTIROSO. No me digas esto y luego me utilices el activate para inicializar el estado

	// Si las casillas aledañas no estaban ocupadas ahora sí lo estarán porque tendrán una casilla vecina que sí lo está
	if (p.i > 0 && !is_occ(Position(p.i - 1, p.j))) {
                config_occ.push_back(Position(p.i - 1, p.j));
                N_occ++;
        }
        if (p.i < config.size() - 1 && !is_occ(Position(p.i + 1, p.j))) {
                config_occ.push_back(Position(p.i + 1, p.j));
                N_occ++;
        }
        if (p.j > 0 && !is_occ(Position(p.i, p.j - 1))) {
                config_occ.push_back(Position(p.i, p.j - 1));
                N_occ++;
        }
        if (p.j < config[0].size() - 1 && !is_occ(Position(p.i, p.j + 1))) {
                config_occ.push_back(Position(p.i, p.j + 1));
                N_occ++;
        }

	auto it = find(config_occ.cbegin(), config_occ.cend(), p);
	if (it == config_occ.cend()) {
		config_occ.push_back(p);
		added = true;
		N_occ++;
	}

	W += (added ? 1 : 1 - get_n(p)*landa/2/d);
	if (p.i > 0 && config[p.i-1][p.j] == 0) W += landa/2/d;
        if (p.i < config.size() - 1 && config[p.i+1][p.j] == 0) W += landa/2/d;
        if (p.j > 0 && config[p.i][p.j-1] == 0) W += landa/2/d;
        if (p.j < config.size() - 1 && config[p.i][p.j+1] == 0) W += landa/2/d;

	config[p.i][p.j] = 1;
}

void CPSim::deactivate(const Position& p) {
	if (config[p.i][p.j] == 0) return; // Si la casilla ya estaba desactivada vamos a cubrirnos las espaldas

	// Si la casilla se propuso para desactivarse tanto ella como sus vecinas estarán en config_occ

	config[p.i][p.j] = 0;

	W += get_n(p)*landa/2/d - 1;
	if (p.i > 0 && config[p.i-1][p.j] == 0) W -= landa/2/d;
	if (p.i < config.size() - 1 && config[p.i+1][p.j] == 0) W -= landa/2/d;
	if (p.j > 0 && config[p.i][p.j-1] == 0) W -= landa/2/d;
	if (p.j < config.size() - 1 && config[p.i][p.j+1] == 0) W -= landa/2/d;

	if (!is_occ(p)) {
		erase_occ_pos(p);
	}
	if (p.i > 0 && !is_occ(Position(p.i-1, p.j))) {
		erase_occ_pos(Position(p.i-1, p.j));
	}
	if (p.i < config.size() - 1 && !is_occ(Position(p.i+1, p.j))) {
		erase_occ_pos(Position(p.i+1, p.j));
	}
	if (p.j > 0 && !is_occ(Position(p.i, p.j-1))) {
		erase_occ_pos(Position(p.i, p.j-1));
	}
	if (p.j < config.size() - 1 && !is_occ(Position(p.i, p.j+1))) {
		erase_occ_pos(Position(p.i, p.j+1));
	}
}

void CPSim::update() {
	if (N_occ == 0) return;

	double p, r;
	t += -log(ran_u(gen))/W;

	auto it = config_occ.cbegin();
	p = ran_u(gen)*W;
	r = (config[it->i][it->j] == 0 ? get_n(*it)*landa/2/d : 1);
	while (r < p && ++it != config_occ.cend()) {
		r += (config[it->i][it->j] == 0 ? get_n(*it)*landa/2/d : 1);
	}

	if (it != config_occ.cend()) {
		if (config[it->i][it->j] == 0) fast_activate(*it);
		else deactivate(*it);
	}
}

void simula_landas(int N, double tmax, int start, int stop, const string& filename) {
	double t_save, freq_t = 1;
	double landa;

	int i, j, aux, N_rho;
	double rho;
	ofstream fout(filename.c_str());

	CPSim simulacion;

	for (i = start; i < stop; i++) {
		landa = 1.6488 + pow(1.13, -i);
		{
			osyncstream out{cout};
			out << "Lanzando simulaciones para landa = " << landa << '\n';

		}
		fout << landa;

		j = aux = 0;
		while (j < 10) {
			simulacion = CPSim(N, landa, 0, SEED*(i+10)*(aux+1));
			simulacion.set_random_config(.5);
			rho = t_save = .0; N_rho = 0;

			while (simulacion.get_t() < tmax && simulacion.get_N_occ() != 0) {
				if (t_save < simulacion.get_t()) {
					if (simulacion.get_t() > 200) {
						rho += simulacion.get_rho();
						N_rho++;
					}
					t_save += freq_t;
				}
				simulacion.update();
			}

			if (simulacion.get_N_occ() != 0) {
				j++;
				fout << '\t' << rho/N_rho;
			}
			aux++;
		}
		fout << endl;
	}

	fout.close();
}

int main() {
	int N = 100;
	double tmax = 500;

	jthread t1(simula_landas, N, tmax, -5, 4, "2D_rho_vs_landa_grandes.txt");
	jthread t2(simula_landas, N, tmax, 4, 13, "2D_rho_vs_landa_medianos.txt");
	jthread t3(simula_landas, N, tmax, 13, 23, "2D_rho_vs_landa_pequenyos.txt");
}
