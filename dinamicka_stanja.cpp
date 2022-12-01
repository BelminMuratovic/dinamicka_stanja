#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

class Motor
{
private:
    double Pn_ = 0;        // nominalna snaga
    double Un_ = 0;        // nominalni napon
    double In_ = 0;        // nominalna struja
    double n_n_ = 0;       // nominalna brzina obrtanja
    double n_s_ = 0;       // sinhrona brzina obrtanja
    double Mn_ = 0;        // nominalni obrtni moment
    double IpIn_ = 0;      // odnos polazne i nominalne struje
    double MpMn_ = 0;      // odnos polaznog i nominalnog momenta
    double MprMn_ = 0;     // odnos prekretnog i nominalnog momenta
    double Jrm_ = 0;       // moment inercije rotora
    double cos_fi_m_ = 0;  // nominalni faktor snage
    double cos_fi_mp_ = 0; // faktor snage pri pokretanju
    double eta_ = 0;       // stepen korisnog dejstva
    double R1_ = 0;        // otpor po fazi statora
    double m1_ = 0;        // masa bakra statora
    double m2_vk_ = 0;     // masa bakra rotora (vanjski kavez)
    double c1_ = 0;        // specificna toplota bakra statora
    double c2_ = 0;        // specificna toplota bakra rotora
    double k1_ = 0;        // koeficijent odvodjenja toplote sa statora
    double k2_ = 0;        // koeficijent odvodjenja toplote sa rotora
    double k3_ = 0;        // koeficijent raspodjele gubitaka izmedju kaveza
    double t_ = 0;         // temperatura radne okoline
    double s_n_ = 0;       // nominalno klizanje

public:
    void klizanje() { s_n_ = (n_s_ - n_n_) / n_s_; }

    double vrijeme_pokretenja_EMP(double Jt, double M)
    {
        return (((Jt + Jrm_) * n_s_) / 9.55) * (1 / M) * (1 - s_n_);
    }

    double gubici_u_bakru_statora(double Jt, double I, double M)
    {
        return (3 * R1_ * ((Jt + Jrm_) * n_s_) / 9.55) * (pow(I, 2) / M) * (1 - s_n_);
    }

    double nadtemperatura_statorskog_namotaja(double Acu1)
    {
        return (Acu1 * k1_) / (m1_ * c1_);
    }

    double temperatura_namotaja_statora(double t)
    {
        return t + t_;
    }

    double gubici_u_bakru_rotora(double Jt, double Mm, double M)
    {
        return ((Jt + Jrm_) * pow(n_s_, 2) / 91.2) * (Mm / (M * 2)) * (1 - pow(s_n_, 2));
    }

    double nadtemperatura_bakra_vanjskog_kaveza(double Acu2)
    {
        return (Acu2 * k2_ * k3_) / (m2_vk_ * c2_);
    }

    double temperatura_vanjskog_kaveza(double t)
    {
        return t + t_;
    }

    friend std::istream &operator>>(std::istream &in, Motor &motor);
};

void unos_podataka(Motor &motor, double &teret)
{
    std::cout << "Unesite specifikacije:" << std::endl;
    std::cout << std::endl;
    std::cin >> motor;
    std::cout << "Moment inercije tereta: ";
    std::cin >> teret;
}

std::istream &operator>>(std::istream &in, Motor &motor)
{
    std::cout << "Nominalna snaga (u kW): ";
    in >> motor.Pn_;
    std::cout << "Nominalni napon (u V): ";
    in >> motor.Un_;
    std::cout << "Nominalna struja (u A): ";
    in >> motor.In_;
    std::cout << "Nominalna brzina obrtanja: ";
    in >> motor.n_n_;
    std::cout << "Sinhrona brzina obrtanja: ";
    in >> motor.n_s_;
    std::cout << "Nominalni obrtni moment (u Nm): ";
    in >> motor.Mn_;
    std::cout << "Odnos polazne i nominalne struje: ";
    in >> motor.IpIn_;
    std::cout << "Odnos polaznog i nominalnog momenta: ";
    in >> motor.MpMn_;
    std::cout << "Odnos prekretnog i nominalnog momenta: ";
    in >> motor.MprMn_;
    std::cout << "Moment inercije rotora: ";
    in >> motor.Jrm_;
    std::cout << "Nominalni faktor snage: ";
    in >> motor.cos_fi_m_;
    std::cout << "Faktor snage pri pokretanju: ";
    in >> motor.cos_fi_mp_;
    std::cout << "Stepen korisnog dejstva: ";
    in >> motor.eta_;
    std::cout << "Otpor po fazi statora: ";
    in >> motor.R1_;
    std::cout << "Masa bakra statora: ";
    in >> motor.m1_;
    std::cout << "Masa bakra rotora (vanjski kavez): ";
    in >> motor.m2_vk_;
    std::cout << "Specificna toplota bakra statora: ";
    in >> motor.c1_;
    std::cout << "Specificna toplota bakra rotora: ";
    in >> motor.c2_;
    std::cout << "Koeficijent odvodjenja toplote sa statora: ";
    in >> motor.k1_;
    std::cout << "Koeficijent odvodjenja toplote sa rotora: ";
    in >> motor.k2_;
    std::cout << "Koeficijent raspodjele gubitaka izmedju kaveza: ";
    in >> motor.k3_;
    std::cout << "Temperatura radne okoline (u °C): ";
    in >> motor.t_;

    return in;
}

void unos_momenata_motora(std::vector<std::pair<double, double>> &momenti)
{
    std::ifstream in{"moment_motora.txt"};
    if (in.is_open())
    {
        while (!in.fail())
        {
            double __s, __M;
            in >> __s >> __M;
            momenti.push_back({__s, __M});
        }
        in.close();
    }
    else
    {
        std::cout << "Fajl se ne nalazi u odgovarajucem folderu" << std::endl;
    }
}

void unos_momenata_tereta(std::vector<std::pair<double, double>> &momenti)
{
    std::ifstream in{"moment_tereta.txt"};
    if (in.is_open())
    {
        while (!in.fail())
        {
            double __s, __M;
            in >> __s >> __M;
            momenti.push_back({__s, __M});
        }
        in.close();
    }
    else
    {
        std::cout << "Fajl se ne nalazi u odgovarajucem folderu" << std::endl;
    }
}

void unos_struje(std::vector<std::pair<double, double>> &struja)
{
    std::ifstream in{"struja.txt"};
    if (in.is_open())
    {
        while (!in.fail())
        {
            double __s, __I;
            in >> __s >> __I;
            struja.push_back({__s, __I});
        }
        in.close();
    }
    else
    {
        std::cout << "Fajl se ne nalazi u odgovarajucem folderu" << std::endl;
    }
}

double povrsina_ispod_krive(const std::vector<std::pair<double, double>> &baza_podataka)
{
    double povrsina = 0;
    double povrsina_odjeljka = 0;
    double gornja_granica = 0;

    if (baza_podataka.size() > 1)
    {
        for (int i = 0; i < baza_podataka.size() - 1; i++)
        {
            gornja_granica = (baza_podataka[i].second + baza_podataka[i + 1].second) / 2;
            povrsina_odjeljka = gornja_granica * (baza_podataka[i].first - baza_podataka[i + 1].first);
            povrsina += povrsina_odjeljka;
        }
    }
    else
    {
        std::cout << "Nedovoljno podataka za proracun" << std::endl;
    }
    return povrsina;
}

double razlika_momenata(const std::vector<std::pair<double, double>> &moment_motora, const std::vector<std::pair<double, double>> &moment_tereta)
{
    return povrsina_ispod_krive(moment_motora) - povrsina_ispod_krive(moment_tereta);
}

void ispis_rezultata_proracuna(double vrijeme_pokretanja, double gubici_u_bakru_statora, double nadtemperatura_statorskog_namotaja, double temperatura_namotaja_statora, double gubici_u_bakru_rotora, double nadtemperatura_bakra_vanjskog_kaveza, double temperatura_vanjskog_kaveza)
{
    std::cout << "\n"
              << std::endl;
    std::cout << "REZULTATI PRORACUNA" << std::endl;
    std::cout << std::endl;
    std::cout << "Vrijeme pokretanja EMP                                 t = " << vrijeme_pokretanja << "(s)" << std::endl;
    std::cout << "Gubici u bakru statora                                 Acu1 = " << gubici_u_bakru_statora / 1000 << "(kJ)" << std::endl;
    std::cout << "Nadtemperatura statorskog namotaja                     Δθ1 = " << nadtemperatura_statorskog_namotaja << "(°C)" << std::endl;
    std::cout << "Temperatura statorskog namotaja                        θ1 = " << temperatura_namotaja_statora << "(°C)" << std::endl;
    std::cout << "Gubici u bakru rotora                                  Acu2 = " << gubici_u_bakru_rotora / 1000 << "(kJ)" << std::endl;
    std::cout << "Nadtemperatura rotorskog namotaja (vanjski kavez)      Δθ2vk = " << nadtemperatura_bakra_vanjskog_kaveza << "(°C)" << std::endl;
    std::cout << "Temperatura rotorskog namotaja (vanjski kavez)         θ2vk = " << temperatura_vanjskog_kaveza << "(°C)" << std::endl;
    std::cout << std::endl;
}

int main()
{
    Motor motor;
    double ukupna_struja = 0;
    double moment_motora_Mm = 0;
    double moment_inercije_tereta = 0;
    double razlika_momenata_Mm_Mt = 0;
    std::vector<std::pair<double, double>> struja;
    std::vector<std::pair<double, double>> moment_motora;
    std::vector<std::pair<double, double>> moment_tereta;

    double vrijeme_pokretanja = 0;
    double gubici_u_bakru_rotora = 0;
    double gubici_u_bakru_statora = 0;
    double temperatura_vanjskog_kaveza = 0;
    double temperatura_namotaja_statora = 0;
    double nadtemperatura_statorskog_namotaja = 0;
    double nadtemperatura_bakra_vanjskog_kaveza = 0;

    unos_podataka(motor, moment_inercije_tereta);
    motor.klizanje();
    unos_momenata_motora(moment_motora);
    unos_momenata_tereta(moment_tereta);
    unos_struje(struja);

    razlika_momenata_Mm_Mt = razlika_momenata(moment_motora, moment_tereta);
    moment_motora_Mm = povrsina_ispod_krive(moment_motora);
    ukupna_struja = povrsina_ispod_krive(struja);

    vrijeme_pokretanja = motor.vrijeme_pokretenja_EMP(moment_inercije_tereta, razlika_momenata_Mm_Mt);
    gubici_u_bakru_statora = motor.gubici_u_bakru_statora(moment_inercije_tereta, ukupna_struja, razlika_momenata_Mm_Mt);
    nadtemperatura_statorskog_namotaja = motor.nadtemperatura_statorskog_namotaja(gubici_u_bakru_statora);
    temperatura_namotaja_statora = motor.temperatura_namotaja_statora(nadtemperatura_statorskog_namotaja);
    gubici_u_bakru_rotora = motor.gubici_u_bakru_rotora(moment_inercije_tereta, moment_motora_Mm, razlika_momenata_Mm_Mt);
    nadtemperatura_bakra_vanjskog_kaveza = motor.nadtemperatura_bakra_vanjskog_kaveza(gubici_u_bakru_rotora);
    temperatura_vanjskog_kaveza = motor.temperatura_vanjskog_kaveza(nadtemperatura_bakra_vanjskog_kaveza);

    ispis_rezultata_proracuna(vrijeme_pokretanja, gubici_u_bakru_statora, nadtemperatura_statorskog_namotaja, temperatura_namotaja_statora, gubici_u_bakru_rotora, nadtemperatura_bakra_vanjskog_kaveza, temperatura_vanjskog_kaveza);

    return 0;
}
