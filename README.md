# Reproduksi Simulasi N-Body Paralel

Dokumen ini menjelaskan prosedur reproduksi seluruh data mentah hasil ukur pada
direktori `results/` untuk studi kasus simulasi gravitasi N-body. Simulasi
dijalankan sebagai solusi paralel berbasis OpenMPI di dalam container Apptainer
pada cluster HPC dua node, menggunakan `mpirun` multi-node karena pengendali Slurm
tidak tersedia. Skrip Slurm tersedia untuk portabilitas.

## 1. Formulasi Masalah dan Algoritma Serial

Percepatan setiap benda dihitung dari interaksi dengan seluruh benda lain,

$$
\mathbf{a}_i = G \sum_{j \ne i} m_j
\frac{\mathbf{r}_j-\mathbf{r}_i}{\left(\lVert\mathbf{r}_j-\mathbf{r}_i\rVert^2+\epsilon^2\right)^{3/2}},
$$

dengan satuan AU, AU/hari, dan massa Matahari, serta \(G=k^2\),
\(k=0.01720209895\), dan \(\epsilon^2=10^{-12}\). Integrasi waktu menggunakan
Velocity Verlet bentuk kick-drift-kick. Setiap langkah melakukan dua evaluasi
gaya, sehingga kompleksitas waktu algoritma serial adalah \(O(N^2)\) per langkah
dengan kebutuhan memori \(O(N)\).

Versi serial pada `src/nbody_serial.c` berperan sebagai baseline. Kernel fisika
yang sama digunakan bersama oleh seluruh varian melalui `src/physics.c`.

## 2. Desain Paralel

Strategi dekomposisi yang digunakan adalah dekomposisi data pada indeks target
secara kontigu (dekomposisi baris satu dimensi). Jumlah target per rank berbeda
paling banyak satu, sehingga penanganan load balance bersifat statis dan seimbang.
Setiap rank menyimpan keadaan global, tetapi hanya menghitung gaya untuk target
lokal.

Pola komunikasi seluruhnya menggunakan operasi kolektif. Posisi lokal
dikumpulkan dengan `MPI_Allgatherv` setelah pergeseran posisi, kecepatan lokal
dikumpulkan dengan `MPI_Allgatherv` kedua setelah pembaruan setengah kecepatan
berikutnya, dan energi kinetik serta potensial dijumlahkan dengan
`MPI_Allreduce`. Nilai waktu komputasi dan komunikasi direduksi dengan `MPI_Reduce`
menggunakan operasi maksimum. Tidak terdapat pasangan `MPI_Send`/`MPI_Recv`, sehingga
seluruh rank mengikuti urutan kolektif yang sama pada setiap langkah.

Kerja evaluasi gaya per rank berorde \(O(N^2/P)\), sedangkan setiap operasi
kolektif menerima keadaan global berukuran \(O(N)\). Rasio
komputasi:komunikasi karena itu berorde \(N/P\); nilai aktualnya diamati melalui
kolom `t_comp` dan `t_comm` pada data hasil ukur.

Varian MPI dibangun dalam dua bentuk dari satu berkas sumber
`src/nbody_mpi.c`: `nbody_mpi` (MPI murni) dan `nbody_mpi_omp` (MPI+OpenMP,
yaitu bonus hibrida). Varian OpenMP membagi perulangan target lokal dengan penjadwalan
statis.

## 3. Kontainerisasi dan Model Eksekusi

`container/nbody.def` merupakan resep Apptainer. Resep ini membangun image berbasis
Ubuntu 22.04, memasang OpenMPI
4.1.2 ke `/usr/local`, menyalin `src/` ke `/opt/nbody/src`, lalu membangun ketiga
varian di `/opt/nbody`.

Model eksekusi yang digunakan adalah hibrida: `mpirun` pada host meluncurkan rank
ke setiap node, dan setiap rank menjalankan biner MPI di dalam container melalui
`apptainer exec`. Direktori `/home` dipetakan dengan opsi `--bind` agar jalur data
dan hasil identik di host dan container.

Kompatibilitas versi MPI disyaratkan: OpenMPI host, instalasi ruang pengguna pada
`$HOME/opt/ompi412`, dan OpenMPI di dalam image harus sama-sama 4.1.2.
Ketidakcocokan versi dapat menggagalkan inisiasi proses atau komunikasi antarrank.

## 4. Prasyarat Lingkungan

- Linux dengan `cc`, `mpicc`, `make`, dan OpenMPI 4.1.2.
- `apptainer` dengan dukungan `--fakeroot` serta akses jaringan untuk mengambil
  basis Ubuntu dan kode sumber OpenMPI.
- Cluster dua node, `compute01` sebagai node entri dan `compute02` yang hanya
  dapat dijangkau dari `compute01`, tanpa sistem berkas bersama.
- Python 3 dengan `matplotlib` untuk skrip plotting.

## 5. Struktur Direktori

- `src/` berisi kode sumber C, `make_data.py` untuk membangun ulang data inti, dan
  `Makefile`.
- `container/` berisi `nbody.def` dan image `nbody.sif`.
- `scripts/` berisi `hostfile`, `slurm_nbody.sh`, dan `plot.py` untuk skrip
  plotting.
- `data/` berisi data masukan: keadaan awal tata surya inti, katalog asteroid,
  serta data penskalaan kuat dan lemah.
- `results/` berisi data mentah hasil ukur.
- `fig/` berisi keluaran skrip plotting.

## 6. Kompilasi, Pembangunan Image, dan Sinkronisasi

Seluruh perintah dijalankan dari root directory proyek pada `compute01`, kecuali
dinyatakan lain. Kompilasi menghasilkan biner di `src/bin/`:

```bash
make -C src
```

Image dibangun dari direktori `container/` karena sumber `%files` dihitung
relatif terhadap direktori kerja:

```bash
(cd container && apptainer build --fakeroot nbody.sif nbody.def)
```

Kedua node tidak berbagi sistem berkas, sehingga direktori proyek harus berada pada
jalur absolut yang sama di `compute01` dan `compute02`. Berkas `scripts/hostfile`
mendefinisikan delapan slot pada masing-masing node:

```text
compute01 slots=8
compute02 slots=8
```

Sinkronisasi dan pemeriksaan checksum image:

```bash
ssh compute02 "mkdir -p \"$PWD\""
rsync -avz --delete --exclude .git "$PWD/" "compute02:$PWD/"
md5sum container/nbody.sif
ssh compute02 "md5sum \"$PWD/container/nbody.sif\""
```

Kedua nilai `md5sum` harus identik. Seluruh perintah `mpirun` berikut memerlukan
MPI ruang pengguna pada `PATH` dan `LD_LIBRARY_PATH`, serta jumlah thread OpenMP.
Jalankan pada shell yang sama:

```bash
export PATH="$HOME/.local/bin:$PATH"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
```

## 7. Data Mentah Hasil Ukur

Setiap subbagian menyatakan isi data, konfigurasi eksperimen, dan berkas
pendukungnya pada `results/`.

### 7.1 Waktu Eksekusi Penskalaan Kuat (Strong Scaling)

Penskalaan kuat mempertahankan ukuran masalah \(N=8000\) dan menaikkan jumlah
proses \(P=1,2,4,8,16\). Setiap titik diulang 5 kali dengan 1000 langkah dan
\(\Delta t=0.5\) hari. Data memuat waktu total yang diukur pada host, nilai
maksimum waktu komputasi dan komunikasi antar rank, serta energi akhir setiap ulangan.
Berkas pendukung: `results/strong.csv` dengan kolom
`N,P,trial,t_total,t_comp,t_comm,t_sum,energy`; kolom `trial` menyatakan nomor
ulangan.

```bash
apptainer exec --containall --bind /home:/home container/nbody.sif \
  /opt/nbody/nbody_serial "$PWD/data/strong8000.csv" 1000 0.5 \
  "$PWD/results/ref_snap.csv" "$PWD/results/ref_energy.txt"
echo "N,P,trial,t_total,t_comp,t_comm,t_sum,energy" > results/strong.csv
for P in 1 2 4 8 16; do
  for t in 1 2 3 4 5; do
    start=$(date +%s.%N)
    out=$(mpirun --hostfile scripts/hostfile -np "$P" -x PATH -x LD_LIBRARY_PATH -x OMP_NUM_THREADS --bind-to core \
      apptainer exec --bind /home:/home container/nbody.sif /opt/nbody/nbody_mpi_omp \
      data/strong8000.csv 1000 0.5 results/snap.csv results/energy.txt)
    end=$(date +%s.%N)
    total=$(awk "BEGIN {printf \"%.6f\", $end - $start}")
    tcomp=$(echo "$out" | sed -n 's/.*t_comp=\([^ ]*\) .*/\1/p')
    tcomm=$(echo "$out" | sed -n 's/.*t_comm=\([^ ]*\) .*/\1/p')
    tsum=$(echo "$out" | sed -n 's/.*t_sum=\([^ ]*\).*/\1/p')
    echo "8000,$P,$t,$total,$tcomp,$tcomm,$tsum,$(cat results/energy.txt)" >> results/strong.csv
    python3 -c "import sys; ref=float(open('results/ref_energy.txt').read().split()[0]); new=float(open('results/energy.txt').read().split()[0]); d=abs(ref-new)/max(abs(ref),abs(new),1e-300); print('P=$P ulangan=$t reldiff=%.3e' % d); sys.exit(0 if d <= 1e-12 else 1)" || exit 1
  done
done
```

Durasi ±1,5 jam. Rata-rata dan simpangan baku waktu dari seluruh ulangan dihitung
dari berkas ini, demikian pula speedup, efisiensi, dan fraksi serial Karp-Flatt
melalui `scripts/plot.py`.

### 7.2 Waktu Eksekusi Penskalaan Lemah (Weak Scaling)

Penskalaan lemah mempertahankan beban per proses dengan \(N=1000P\), sehingga
data yang digunakan adalah `data/weak1000.csv` hingga `data/weak16000.csv`.
Setiap titik diulang 3 kali dengan 1000 langkah dan \(\Delta t=0.5\) hari.
Referensi serial dihitung untuk setiap data karena setiap titik memiliki ukuran masalah
berbeda. Berkas pendukung: `results/weak.csv` dengan kolom yang sama seperti
penskalaan kuat.

```bash
for P in 1 2 4 8 16; do
  N=$((1000*P))
  apptainer exec --containall --bind /home:/home container/nbody.sif \
    /opt/nbody/nbody_serial "$PWD/data/weak$N.csv" 1000 0.5 \
    "$PWD/results/ref_snap_$N.csv" "$PWD/results/ref_energy_$N.txt"
done
echo "N,P,trial,t_total,t_comp,t_comm,t_sum,energy" > results/weak.csv
for P in 1 2 4 8 16; do
  N=$((1000*P))
  for t in 1 2 3; do
    start=$(date +%s.%N)
    out=$(mpirun --hostfile scripts/hostfile -np "$P" -x PATH -x LD_LIBRARY_PATH -x OMP_NUM_THREADS --bind-to core \
      apptainer exec --bind /home:/home container/nbody.sif /opt/nbody/nbody_mpi_omp \
      data/weak$N.csv 1000 0.5 results/snap.csv results/energy.txt)
    end=$(date +%s.%N)
    total=$(awk "BEGIN {printf \"%.6f\", $end - $start}")
    tcomp=$(echo "$out" | sed -n 's/.*t_comp=\([^ ]*\) .*/\1/p')
    tcomm=$(echo "$out" | sed -n 's/.*t_comm=\([^ ]*\) .*/\1/p')
    tsum=$(echo "$out" | sed -n 's/.*t_sum=\([^ ]*\).*/\1/p')
    echo "$N,$P,$t,$total,$tcomp,$tcomm,$tsum,$(cat results/energy.txt)" >> results/weak.csv
    python3 -c "import sys; ref=float(open('results/ref_energy_$N.txt').read().split()[0]); new=float(open('results/energy.txt').read().split()[0]); d=abs(ref-new)/max(abs(ref),abs(new),1e-300); print('P=$P ulangan=$t reldiff=%.3e' % d); sys.exit(0 if d <= 1e-12 else 1)" || exit 1
  done
done
```

Durasi ±1 jam. Rata-rata dan simpangan baku waktu dari seluruh ulangan dihitung
dari berkas ini.

### 7.3 Keadaan Akhir Simulasi

Data berisi posisi dan kecepatan setiap benda setelah 1000 langkah pada
\(N=8000\), \(\Delta t=0.5\) hari. Hasil bit-identik untuk \(P\) berapa pun;
\(P=16\) digunakan karena paling cepat (±40 detik). Nama berkas keluaran harus
persis seperti yang tertulis karena dibaca skrip plotting. Berkas pendukung:
`results/snap_final_n8000.csv`.

```bash
mpirun --hostfile scripts/hostfile -np 16 -x PATH -x LD_LIBRARY_PATH -x OMP_NUM_THREADS --bind-to core \
  apptainer exec --bind /home:/home container/nbody.sif /opt/nbody/nbody_mpi_omp \
  data/strong8000.csv 1000 0.5 "$PWD/results/snap_final_n8000.csv" /tmp/energy_final.txt
```

Grafik penskalaan dan visualisasi posisi dihasilkan dengan:

```bash
python3 scripts/plot.py
```

### 7.4 Bukti Validasi Numerik terhadap Baseline Serial

Validasi membandingkan energi setiap ulangan pada penskalaan kuat dan lemah
terhadap referensi serial, dengan toleransi relatif \(10^{-12}\). Berkas
pendukung: `results/validation.log`. Prasyarat: `results/ref_energy.txt`
(subbagian 7.1) dan `results/ref_energy_$N.txt` (subbagian 7.2), sehingga tidak
diperlukan eksekusi serial tambahan. Kedua jenis berkas referensi tersedia,
sehingga validasi dapat dijalankan tanpa menunggu sapuan penskalaan.

```bash
{
echo "== referensi serial vs energi ulangan tercatat =="
for N in 1000 2000 4000 8000 16000; do
  python3 - "$N" <<'PY' || exit 1
import csv, sys
n = sys.argv[1]
worst, count = 0.0, 0
for path in ["results/weak.csv"] + (["results/strong.csv"] if n == "8000" else []):
    ref = float(open("results/ref_energy.txt" if path.endswith("strong.csv")
                     else "results/ref_energy_%s.txt" % n).read().split()[0])
    for row in csv.DictReader(open(path)):
        if int(row["N"]) != int(n):
            continue
        energy = float(row["energy"])
        diff = abs(ref - energy) / max(abs(ref), abs(energy), 1e-300)
        worst = max(worst, diff)
        count += 1
print("N=%s ulangan=%d worst_reldiff=%.3e" % (n, count, worst))
sys.exit(0 if worst <= 1e-12 and count > 0 else 1)
PY
done
echo "ALL VALIDATED"
} 2>&1 | tee results/validation.log
```

Validasi dinyatakan lulus apabila setiap baris `worst_reldiff` tidak melebihi
\(10^{-12}\) dan berkas diakhiri `ALL VALIDATED`.

### 7.5 Waktu Baseline Serial

Data berisi waktu eksekusi baseline serial pada \(N=8000\), 1000 langkah,
\(\Delta t=0.5\) hari, dalam dua lingkungan: bare-metal dan container. Kedua
nilai digunakan sebagai pembanding untuk analisis speedup varian paralel terhadap
varian serial. Berkas pendukung: `results/serial_baseline.csv` dengan kolom
`N,steps,binary,env,t_total,energy`.

```bash
start=$(date +%s.%N)
src/bin/nbody_serial data/strong8000.csv 1000 0.5 results/snap_serial_bare.csv results/energy_serial_bare.txt
end=$(date +%s.%N)
bare=$(awk "BEGIN {printf \"%.3f\", $end - $start}")
start=$(date +%s.%N)
apptainer exec --containall --bind /home:/home container/nbody.sif \
  /opt/nbody/nbody_serial "$PWD/data/strong8000.csv" 1000 0.5 \
  "$PWD/results/snap_serial_container.csv" "$PWD/results/energy_serial_container.txt"
end=$(date +%s.%N)
container=$(awk "BEGIN {printf \"%.3f\", $end - $start}")
echo "N,steps,binary,env,t_total,energy" > results/serial_baseline.csv
printf "8000,1000,bare:src/bin/nbody_serial,bare-metal@%s,%s,%s\n" \
  "$(hostname)" "$bare" "$(cat results/energy_serial_bare.txt)" >> results/serial_baseline.csv
printf "8000,1000,container:/opt/nbody/nbody_serial,apptainer+nbody.sif@%s,%s,%s\n" \
  "$(hostname)" "$container" "$(cat results/energy_serial_container.txt)" >> results/serial_baseline.csv
```

Durasi ±8 menit per eksekusi. Kolom `binary` mencatat jalur biner saat pengukuran.

### 7.6 Dampak Kontainerisasi terhadap Kinerja

Data berisi perbandingan berpasangan mode bare-metal dan container pada
konfigurasi tetap \(N=8000\), 200 langkah, \(P=4\), dan 3 ulangan per mode, yang
seluruhnya dijalankan pada `compute01`. Mode bare-metal menjalankan
`src/bin/nbody_mpi_omp` secara langsung, sedangkan mode container menjalankan
biner yang sama di dalam image. Data ini mendukung bonus perbandingan overhead
container terhadap bare-metal. Berkas pendukung: `results/overhead.csv` dengan
kolom `N,steps,P,trial,machine,config,t_total,energy`.

```bash
echo "N,steps,P,trial,machine,config,t_total,energy" > results/overhead.csv
for t in 1 2 3; do
  for mode in bare container; do
    start=$(date +%s.%N)
    if [ "$mode" = bare ]; then
      mpirun -np 4 --bind-to core src/bin/nbody_mpi_omp \
        data/strong8000.csv 200 0.5 results/snap.csv results/energy.txt > /dev/null
    else
      mpirun -np 4 --bind-to core \
        apptainer exec --bind /home:/home container/nbody.sif /opt/nbody/nbody_mpi_omp \
        data/strong8000.csv 200 0.5 results/snap.csv results/energy.txt > /dev/null
    fi
    end=$(date +%s.%N)
    total=$(awk "BEGIN {printf \"%.3f\", $end - $start}")
    echo "8000,200,4,$t,$(hostname),$mode,$total,$(cat results/energy.txt)" >> results/overhead.csv
  done
done
```

Durasi ±5 menit. Pastikan `OMP_NUM_THREADS=1` dari bagian 6 masih berlaku agar
perbandingan kedua mode setara.

### 7.7 Profiling Program Serial dengan gprof

Data berisi distribusi sampel waktu per fungsi pada varian serial, diambil
dengan `gprof`, sebagai bukti bonus profiling. Berkas pendukung:
`results/gprof.txt`. Profil yang dihasilkan menempatkan 99,54% sampel waktu pada
`physics_forces`.

```bash
cc -O3 -march=native -std=c11 -Wall -pg -o /tmp/nbody_prof src/nbody_serial.c src/loader.c src/physics.c -lm
/tmp/nbody_prof data/strong8000.csv 50 0.5 /tmp/p_snap.csv /tmp/p_energy.txt
gprof /tmp/nbody_prof gmon.out > results/gprof.txt
rm -f gmon.out
```

Beban 50 langkah pada \(N=8000\) digunakan agar sampel waktu cukup terkumpul;
menjalankan beban terlalu kecil menghasilkan profil kosong.

### 7.8 Percobaan Pencacah Perangkat Keras

Data berisi transkrip perintah `perf stat` beserta pesan kegagalannya pada kernel
cluster ini, sebagai catatan keterbatasan pemrofilan perangkat keras.
Berkas pendukung: `results/perf_serial.txt`.

```bash
{ echo "\$ perf stat -e cycles,instructions,cache-misses src/bin/nbody_serial data/solar_system.csv 50 0.5 /tmp/p.csv /tmp/p.txt"
  perf stat -e cycles,instructions,cache-misses src/bin/nbody_serial data/solar_system.csv 50 0.5 /tmp/p.csv /tmp/p.txt && code=0 || code=$?
  echo "perf_exit=$code"
} > results/perf_serial.txt 2>&1
```

Apabila `perf` tersedia pada mesin lain, berkas akan berisi angka pencacah nyata,
bukan pesan kegagalan.

## 8. Eksekusi di Cluster melalui Slurm

Pada cluster yang digunakan, pengendali Slurm tidak tersedia sehingga eksekusi
dilakukan melalui perintah `mpirun` multi-node yang setara, seperti pada bagian 7.
`scripts/slurm_nbody.sh` disediakan untuk cluster yang memiliki pengendali
Slurm, menjalankan satu titik melalui `srun` dan Apptainer. Parameter `DATA`,
`STEPS`, `DT`, `SNAP_OUT`, `ENERGY_OUT`, dan `OMP_NUM_THREADS` diatur melalui
variabel lingkungan, misalnya:

```bash
DATA="$PWD/data/strong8000.csv" STEPS=100 DT=0.5 sbatch scripts/slurm_nbody.sh
```

## 9. Catatan Reproduksibilitas

- Data pada `data/` tersedia agar seluruh prosedur dapat dijalankan tanpa
  membangun ulang data. `src/make_data.py` dapat membangun ulang data inti;
  pembuatan ulang `data/asteroids.csv` memerlukan berkas MPCORB.
- `data/weak8000.csv` dan `data/strong8000.csv` berisi himpunan benda yang sama.
- Seluruh hasil numerik bersifat deterministik; validasi terhadap baseline serial
  menggunakan toleransi relatif energi \(10^{-12}\).
- Waktu eksekusi per rank tidak direkam, sehingga analisis load imbalance dibatasi
  pada nilai maksimum antar rank yang tersedia melalui `t_comp` dan `t_comm`.
