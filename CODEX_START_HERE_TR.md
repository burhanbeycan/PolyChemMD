# PolyChemMD: Codex'te geliştirmeye başlama

Bu paket, bağımsız bir genel polimer hesaplama platformu için başlangıç deposudur. Mevcut hafif çekirdek; kayıt doğrulama, hesaplama kökeni ve açık konvansiyonlu genel Fourier torsiyon analizini kapsar. Polimer oluşturucu, üretim kuvvet alanı, QM/MD motorları ve doğrulanmış malzeme özellikleri henüz tamamlanmış değildir.

## 1. Depoyu GitHub'a yerleştirin

1. GitHub hesabınızda boş bir `PolyChemMD` deposu oluşturun. Çalışma adı daha sonra değiştirilebilir. İlk geliştirme için özel depo kullanılabilir.
2. Paketi açın; `README.md`, `AGENTS.md`, `pyproject.toml`, `src`, `tests` ve `docs` aynı depo kökünde olsun. İç içe ikinci bir proje klasörü oluşturmayın.
3. Dosyaları bu depoya aktarın ve ilk commit'i oluşturun. Paketin paylaşılması, GitHub deposunun oluşturulduğu veya Codex'e bağlandığı anlamına gelmez.
4. Büyük yörüngeleri, motor kurulumlarını ve kişisel erişim bilgilerini Git'e koymayın. Yeniden üretim reçetelerini ve küçük doğrulama örneklerini saklayın.

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
```

Kurulum ortamı ağ erişimine ihtiyaç duyuyorsa kurulum adımında uygun bağımlılıklar sağlanmalıdır. Derleme bağımlılıkları zaten kurulu olan ağsız bir geliştirme ortamında şu kurulum kullanılabilir:

```bash
python -m pip install --no-deps --no-build-isolation -e .
```

`validate` mevcut kayıt şemasını denetler; SMILES'ın moleküler valansını veya polimerin fiziksel doğruluğunu henüz doğrulamaz. `catalogue` planlanan özellikleri de listeler; listedeki adlar hesaplanabilirlik kanıtı değildir. `torsion-describe` analitik örneği işler ve malzeme MD verisi üretmez. `doctor` seçilen makinenin durumunu gösterir; bu kontrol de motorların bilimsel doğrulamasının yerini tutmaz.

Ağsız test koşusunda yeni paket indirmeye bağımlı olmayın. README'deki komutlar ve mevcut kod, gelecekteki sürümlerde bu rehbere göre daha güncel olabilir; sürüm değiştiğinde belgeleri birlikte güncelleyin.

## 5. İşleri küçük ve ölçülebilir görevlere bölün

Önce `docs/TASKS.md` içindeki kimyasal sözleşme ve atom eşleme görevi, ardından doğrulanmış lineer zincir oluşturucu gelir. Sonraki sıra parametre sözleşmesi, gerçek motor bağlantısı ve tek bir uçtan uca doğrulanmış özelliktir. Bütün özellik kataloğunu tek görevde istemeyin.

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
planlananları ayır. TASK-001 için tekrar birimi bağlanma noktaları, atom
eşleme/soy takibi ve uç grup sözleşmesini tamamla. Bu sözleşme hazırsa
TASK-002 kapsamında bağımsız bir lineer tekrar-birimi zincir oluşturucu
uygula. Mevcut koddan farklı bir yaklaşım gerekiyorsa gerekçesini kaydet.

RadonPy veya başka bir polimer paketinin kaynak kodunu kopyalama,
vendor etme veya çalışma zamanı bağımlılığı yapma. Genel kütüphaneler
ve bilimsel yöntemler kaynakları belirtilerek kullanılabilir.

Serbest monomerden polimer bağını oluşturan reaksiyon dönüşümü ile zaten
tanımlanmış tekrar birimlerini bağlamayı ayrı tut. Atom kimlikleri,
bağ derecesi, formal yük, uç gruplar ve stereokimya korunmalı.

Tek, iki ve birkaç tekrar birimi için elle kontrol edilmiş bağlantı ve
atom sayısı örnekleri kullan. Geçersiz valans, tekrar eden atom eşlemesi
ve uyumsuz bağlanma noktaları için anlamlı hatalar ve bağımsız testler yaz.
Opsiyonel RDKit desteği hafif çekirdeğin içe aktarılmasını bozmamalı.

İlgili testleri gerçekten çalıştır ve komutları/sonuçları bildir.
Bu görevde MD/QM çalıştırılmış gibi davranma; koordinat oluşturma veya
fiziksel özellik doğruluğunu kanıtlanmış sayma. Değişiklikleri,
doğrulamayı, sınırlamaları ve sonraki tek görevi raporla.
```

## 7. Her değişiklikten sonra bilimsel inceleme

Kodun çalışması kimyasal doğruluğun ilk kontrolüdür. Bağ yapısı, yük, birim, faz ve istatistiksel örnekleme ayrıca incelenmeli. Yapay analitik örnek, motor testi, sayısal yakınsama ve deneysel karşılaştırma farklı kanıt düzeyleridir.

Makale için katkı iddiasını; kimyasal çevre/mikroyapı duyarlılığı, kuantumdan uzun zincire aktarım ve karşılaştırmayla sınanmış descriptor/özellik yöntemleri üzerine kuracağız. Deponun büyüklüğü veya özellik adlarının sayısı tek başına bilimsel sonuç değildir.
