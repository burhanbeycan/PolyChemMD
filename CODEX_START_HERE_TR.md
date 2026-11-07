# PolyChemMD: Codex'te geliştirmeye başlama

Bu depo, bağımsız bir genel polimer hesaplama platformunun geliştirme çekirdeğidir. Sürüm 0.5; nötr lineer alkanlar için periyodik parametre birleştirme, açık Ewald modeli, bağımsız enerji/kuvvet karşılaştırması ve sınırlı gevşetme–NVT–NPT hazırlığı içerir. Genel kimyasal kapsam, yeterli üretim örneklemesi ve doğrulanmış polimer MD/QM özellikleri henüz tamamlanmış değildir. Kısa hazırlık koşuları denge veya yoğunluk doğrulaması sayılmaz.

## 1. Mevcut depoyla çalışın

Kaynak depo: https://github.com/burhanbeycan/PolyChemMD . Başlangıç paketi bu
depoya aktarılmıştır. Yeni geliştirmeleri mevcut dosyalardan ve güncel ana daldan
başlatın; ikinci bir proje kökü oluşturmayın. Değişiklikleri görev dallarında
inceleyip testleri geçtikten sonra ana dala birleştirin.

Büyük yörüngeleri, motor kurulumlarını ve kişisel erişim bilgilerini Git'e
koymayın. Yeniden üretim reçetelerini ve küçük doğrulama kayıtlarını saklayın.

## 2. Codex'e aynı projeyi açtırın

Kullandığınız arayüzde yerel proje klasörünü açın veya desteklenen GitHub/depo bağlantısı üzerinden bu depoyu seçin. Depo erişimini doğru hesaba ve projeye verin. Ekranların adları arayüze göre değişebilir; gerekli nesneler proje/depo, çalışma dalı ve çalıştırma ortamıdır.

6 Kasım 2026 tarihli resmî Quickstart açıklaması, uzaktan geliştirme için yeni
görevde `Work in > Cloud`, ardından `Select environment > Create environment`
yolunu tarif ediyor. Bu akış hesabınızda görünüyorsa depoyu seçtirin, aşağıdaki
kurulum ve test komutlarıyla ortamı hazırlatın ve tekrar kullanılabilir ortamı
hazır olduğunda yayımlayın. Yerel geliştirmede ise proje klasörünü açıp Codex
seçeneğiyle çalışabilirsiniz. Ortam yapılandırmasının yayımlanması, kod deposunu
kamuya açmak veya bilimsel sonuç yayımlamak anlamına gelmez.

Codex'in AGENTS.md dosyasını ve `docs/SCIENTIFIC_SPEC.md` belgesini okumasını sağlayın. Bu dosyalar; bağımsız uygulama, bilimsel kapsam, gerçek hesaplama ile örnek verinin ayrılması ve kabul ölçütlerini kalıcı hale getirir. Bütün geçmiş sohbeti her göreve yapıştırmak yerine kararları bu dosyalarda güncel tutun.

Resmî açıklamalar:

- Proje talimatları: https://learn.chatgpt.com/docs/agent-configuration/agents-md
- Ayrı çalışma ağaçları: https://learn.chatgpt.com/docs/environments/git-worktrees
- Başlangıç: https://learn.chatgpt.com/docs/quickstart
- Eski Codex bulut ortamı belgeleri: https://learn.chatgpt.com/docs/environments/cloud-environment

Son bağlantı legacy ortamı açıklar; hesabınızda veya kullandığınız arayüzde belirli bir ortamın mevcut olduğunu varsaymayın.

## 3. İki çalıştırma ortamı kullanın

| Ortam | Kullanım | Temel ilke |
| --- | --- | --- |
| Hafif geliştirme | Şema, algoritma, descriptor, CLI, küçük testler ve kod inceleme | MD/QM veya GPU olmadan çalışabilmeli |
| Bilimsel hesaplama | Gerçek polimer hazırlama, LAMMPS ve QM kampanyaları | Motor, parametre, donanım ve örnekleme koşulları kaydedilmeli |

RTX 4090 bulunan makinenizde uygun motor derlemesi ve GPU desteği ayrıca doğrulandıktan sonra ağır hesaplar çalıştırılabilir. Kodlama ortamının kendiliğinden bu GPU'ya veya uygun bilimsel yazılımlara eriştiğini varsaymayacağız. İlk üretim bütçesini kısa zamanlama ve yakınsama pilotlarından sonra belirleyeceğiz.

Temel kütüphane bir OpenAI API anahtarı gerektiren bir uygulama olarak tasarlanmıyor. İleride isteğe bağlı bir LLM ajanı eklenirse onun erişimi ve maliyeti ayrı bir bileşen olur.

## 4. İlk yerel kontrol

Depo kökünde:

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
polychemmd doctor
polychemmd validate examples/pe_record.json
polychemmd catalogue
polychemmd torsion-describe examples/torsion_demo.json --n-grid 720 --temperature-k 300 --output artifacts/torsion_demo.json
polychemmd parameter-check examples/parameter_contract_demo.json --require-complete
polychemmd torsion-import examples/torsion_lammps_fourier.json --source lammps-fourier
```

Kimyasal grafik oluşturmayı denemek için `python -m pip install -e '.[chemistry]'`
ile isteğe bağlı RDKit bağımlılığını kurun; ardından
`polychemmd build-linear examples/linear_pe.json` çalıştırın. Bu komut bir
moleküler grafik üretir; üç boyutlu bir polimer veya simülasyon hücresi üretmez.

Koordinat ve hücre hazırlama komutları ayrıdır:

```bash
polychemmd build-conformer examples/conformer_pe.json --output artifacts/pe-conformer.json
polychemmd pack-cell examples/packed_pe.json --output artifacts/pe-cell.json
polychemmd export-xyz artifacts/pe-cell.json --output artifacts/pe-cell.xyz
```

Bu yapılar başlangıç geometrileridir; üretim kuvvet alanı atanmış veya
dengelenmiş bir polimer fazı değildir. Hücre yoğunluğu seçilen hacimden hesaplanır.

Desteklenen lineer alkan geometrisi için parametre atama ve motor karşılaştırması:

```bash
polychemmd parameterize-alkane artifacts/pe-conformer.json --output artifacts/assigned-pe.json
polychemmd alkane-energy artifacts/pe-conformer.json --output artifacts/pe-energy.json
polychemmd alkane-lammps artifacts/pe-conformer.json --run-dir artifacts/pe-check-plan --plan
```

Son komuttaki `--plan` motoru çalıştırmaz. Gerçek karşılaştırma için LAMMPS
kurulmuş olmalı; `--plan` kaldırıldığında yeni bir koşu dizininde sıfır entegrasyon
adımlı tek-nokta hesabı yürütülür. Bu kontrol yalıtılmış ve sonlu kesimlidir;
periyodik üretim protokolü veya deneysel doğruluk kanıtı değildir.

Periyodik model ve sınırlı hazırlık için ayrı komutlar:

```bash
polychemmd parameterize-cell artifacts/pe-cell.json --output artifacts/assigned-cell.json
polychemmd periodic-lammps artifacts/assigned-cell.json --run-dir artifacts/periodic-plan --plan
polychemmd periodic-prepare artifacts/assigned-cell.json --protocol examples/periodic_protocol.json --run-dir artifacts/preparation-plan --plan
```

Planlar çalıştırma yapmaz. Gerçek motor koşusu için yeni bir dizinde `--plan`
kaldırılır. Hazırlık çıktısı; minimizasyonun durma sebebini, gerçekleşen adımları,
son hücre/kimlik kontrollerini ve örnekleme durumunu ayrı kaydeder. Sıcaklık veya
basınç hedefinin verilmesi, kısa sürede bu hedeflere ulaşıldığını kanıtlamaz.

Kurulum ortamı ağ erişimine ihtiyaç duyuyorsa kurulum adımında uygun bağımlılıklar sağlanmalıdır. Derleme bağımlılıkları zaten kurulu olan ağsız bir geliştirme ortamında şu kurulum kullanılabilir:

```bash
python -m pip install --no-deps --no-build-isolation -e .
```

`validate` mevcut kayıt şemasını denetler; SMILES'ın moleküler valansını veya polimerin fiziksel doğruluğunu henüz doğrulamaz. `catalogue` planlanan özellikleri de listeler; listedeki adlar hesaplanabilirlik kanıtı değildir. `torsion-describe` analitik örneği işler ve malzeme MD verisi üretmez. `doctor` seçilen makinenin durumunu gösterir; bu kontrol de motorların bilimsel doğrulamasının yerini tutmaz.

Ağsız test koşusunda yeni paket indirmeye bağımlı olmayın. README'deki komutlar ve mevcut kod, gelecekteki sürümlerde bu rehbere göre daha güncel olabilir; sürüm değiştiğinde belgeleri birlikte güncelleyin.

## 5. İşleri küçük ve ölçülebilir görevlere bölün

`docs/TASKS.md` güncel durumun kaynağıdır. TASK-001–005, TASK-006A ve sınırlı n-alkan kapsamındaki TASK-006B/C1 uygulanmıştır. Sonraki adım TASK-006C2 durum/hazırlık/örnekleme protokolünün ölçülebilir kaynak bütçesiyle sınanmasıdır. Bağımsız hücreler ve uygun dış kaynaklarla karşılaştırma tam yoğunluk doğrulamasının parçasıdır. Bütün özellik kataloğunu tek görevde istemeyin.

Her görevde şunlar açık olsun:

- Amaç ve değiştirilebilecek kapsam.
- Kimyasal/matematiksel varsayımlar ve birimler.
- Bağımsız doğrulama örnekleri ve kabul ölçütleri.
- Gerçekte çalıştırılacak testler ve motor hesapları.
- Çıktıda kanıtlanması gereken başarı ve açık kalan sınırlar.

Paralel çalışmada her görev ayrı dal/çalışma ağacı ve ayrı dosya sahipliği kullanmalı. Ortak veri şemasını iki görev aynı anda yeniden tasarlamamalı. Ortak sözleşme sabitlendikten sonra oluşturucu, analiz ve motor bağlantıları ayrılabilir.

## 6. Codex'e verilecek ilk görev

```text
PolyChemMD deposunda AGENTS.md, README.md, docs/ARCHITECTURE.md,
docs/SCIENTIFIC_SPEC.md ve docs/TASKS.md belgelerini oku.

Önce mevcut şema ve çalışan çekirdeği incele; uygulanmış özelliklerle
planlananları ayır. TASK-006C2 kapsamında tek ve açıkça tanımlı kimyasal
sistem/faz/durum için hazırlık ve örnekleme protokolünü nitele. İlk pilotun
düşük yoğunluklu dört oligomer ve 200 adımla sınırlı olduğunu dikkate al.
Hücre/zincir büyüklüğü, başlangıç yoğunluğu, uzun menzilli işlemler,
fiziksel süre, bağımsız tohumlar ve ölçülmüş donanım/zaman bütçesini belirle.

RadonPy veya başka bir polimer paketinin kaynak kodunu kopyalama,
vendor etme veya çalışma zamanı bağımlılığı yapma. Genel kütüphaneler
ve bilimsel yöntemler kaynakları belirtilerek kullanılabilir.

Serbest monomerden polimer bağını oluşturan reaksiyon dönüşümü ile zaten
tanımlanmış tekrar birimlerini bağlamayı ayrı tut. Atom kimlikleri,
bağ derecesi, formal yük, uç gruplar ve stereokimya korunmalı.

Mevcut grafik, geometri, atama, periyodik enerji/kuvvet ve protokol testlerini
koru. Kaynak bütçesine uygun yeniden başlatılabilir durumlar ve özellik-bazlı
tanılar ekle. Bağımsız hücreleri korelasyonlu karelerden ayır. Başarısız ve
yetersiz örneklenmiş koşuları görünür tut. Parametre, kesim, süre ve hücre
büyüklüğü duyarlılığını sınamadan; durum/fazı eşleşen dış karşılaştırma
olmadan yoğunluğu doğrulanmış polimer tahmini olarak sunma.

İlgili testleri gerçekten çalıştır ve komutları/sonuçları bildir.
Yalnız gerçekten çalıştırılan motor hesaplarını bildir; üretim MD/QM
ve fiziksel özellik doğruluğunu kanıtlanmış sayma. Değişiklikleri,
doğrulamayı, sınırlamaları ve sonraki tek görevi raporla.
```

## 7. Her değişiklikten sonra bilimsel inceleme

Kodun çalışması kimyasal doğruluğun ilk kontrolüdür. Bağ yapısı, yük, birim, faz ve istatistiksel örnekleme ayrıca incelenmeli. Yapay analitik örnek, motor testi, sayısal yakınsama ve deneysel karşılaştırma farklı kanıt düzeyleridir.

Makale için katkı iddiasını; kimyasal çevre/mikroyapı duyarlılığı, kuantumdan uzun zincire aktarım ve karşılaştırmayla sınanmış descriptor/özellik yöntemleri üzerine kuracağız. Deponun büyüklüğü veya özellik adlarının sayısı tek başına bilimsel sonuç değildir.

## v0.6 geliştirme sınırı

Yerel ikili NPT checkpoint dosyasından devam etme, önceden belirlenen yoğunluk
örnekleme kontrolleri ve bağımsız hücre/kaynak bütçesi planlaması eklendi.
`density-plan examples/density_campaign.json` simülasyon başlatmaz.
Kısa devam koşuları aynı hücrenin parçalarıdır; bağımsız örnek veya dengeye
ulaşmış polimer yoğunluğu değildir. Sonraki görev TASK-006C2B: 512 atom
sınırını ve hedef sıvı yoğunluğu/kesim yarıçapı uyumsuzluğunu giderecek,
kimlik ve parametre denetimlerini koruyan daha büyük periyodik hücre desteği.
