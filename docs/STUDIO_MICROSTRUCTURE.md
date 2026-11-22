# Kontrollü PVA zinciri ve Ti₃C₂Tₓ yüzeyi

Zincir üreticisi belirli sayıda vinil asetat (OAc) komşuluğuna sahip, uçları tanımlı bir oligomer üretir. Yüzey üreticisi her yüzdeki OH ve F sayılarını koruyarak site komşuluklarını değiştirir. Bu iki özellik aynı yapıda kontrol edilebilir. Yapı, SDF dosyası, atomların tekrar birimi kimlikleri ve kullanılan rastgelelik tohumları birlikte kaydedilir.

## PVA dizilişi ve NMR

DP = N ve hedef hidroliz kesri h için asetat sayısı k = ⌊N(1−h)+½⌋ olur. Açık zincirde N−1 komşu tekrar çifti bulunur. D, yan yana bulunan OAc–OAc çiftlerinin sayısıdır. Uçları ayrı kabul edilen bir zincirde aşağıdaki seçenekler vardır.

| Girdi | Uygulama |
|---|---|
| `acetate_pattern: "dyads"`, `acetate_dyads: 1` | Tam olarak bir OAc–OAc komşuluğu içeren diziler arasından eşit olasılıkla örnekler. Geçersiz D için olası tam sayı aralığını verir. |
| `acetate_pattern: "dyads"`, `acetate_dyad_fraction: 0.15` | Bütün N−1 komşuluk arasındaki OAc–OAc kesrini hedefler. Kısa zincir için en yakın mümkün sayıyı seçer ve kesir hatasını verir. Eşit uzaklıkta büyük sayıyı seçer. |
| `acetate_pattern: "explicit"`, `acetate_positions: [3, 6, 9]` | Asetatı girilen 1 tabanlı konumlara yerleştirir. Konum sayısı hidroliz derecesiyle uyuşmalıdır. |
| `acetate_pattern: "random"` veya `"clustered"` | Sabit asetat sayısıyla serbestçe örnekler veya asetatları tek bir kesintisiz blokta toplar. |

**Örnekleme.** k > 0 için r = k−D asetat bloğu bulunur. k asetatı r pozitif uzunluğa ayırıp blokları N−k OH biriminin oluşturduğu N−k+1 aralıktan r tanesine yerleştiririz. Her yönlü dizi tek bir seçim çiftine karşılık gelir. Böylece bir sınıftaki dizi sayısı C(k−1,r−1) × C(N−k+1,r), örnekleme süresi zincir uzunluğu ile sınırlıdır. DP12 ve üç asetat için D = 0, 1, 2 sınıflarında sırasıyla 120, 90 ve 10 dizi vardır. Üç sınıftan eşit örnek alınan çalışma kontrollü bir karşılaştırmadır. Rastgele sentezlenmiş bir popülasyonu temsil etmek için seçilen sentez modeline veya ölçülen dağılıma uygun sınıf ağırlıkları gerekir.

`polymer.sequence_statistics` OH–OH, iki yöndeki OH–OAc ve OAc–OAc sayılarını, kullanılabilir çift sayısını, asetat bloklarının sayısını ve ortalama uzunluğunu saklar. Sıfır komşuluğu olan tek birimde kesir boş bırakılır. NMR piklerinin integrallerini bu tanımla uyumlu olacak şekilde normalize etmek gerekir. Kısmen hidroliz edilmiş PVA/VAc için deneysel dayanak [Moritani ve Fujiwara, Macromolecules](https://doi.org/10.1021/ma60057a007) çalışmasıdır. Ölçülen dizi kesirleri numuneyi topluca tanımlar; tek oligomer aynı kesirlerle uyumlu birçok mikroyapıdan biridir. DP değeri seçilen sonlu oligomerin uzunluğudur. Arayüzün PVA topluluğu hesabı (`pva_ensemble`), seçilen DP aralığında deneysel Mₙ ve dispersiteyle uyumlu, maksimum entropili bir sayı dağılımından açık dizi reçeteleri örnekler. Üretilen her atomistik kutu kendi zincir şablonunu kullanır. Heterojen bir numune için hesaplar arasında topluluk ağırlıklarını ve zincirler arası değişkenliği değerlendirin.

## Yüzey komşuluğu

Her yüz, merkezi Ti atomları üzerindeki üçgensel sonlanma sitelerinden oluşur. İki yönde periyodik bir sistemde her sitenin altı komşusu vardır. İleri üç bağ gerçek uzaydaki her komşu çiftini bir kez sayar. Küçük tekrarlanan hücrelerde aynı site kimlikleri birden fazla farklı periyodik bağa girer; görüntüler ve aynı siteye dönen bağlar bu sebeple korunur.

Nₛ site, E = 3Nₛ bağ, OH kesri f ve OH–OH bağ sayısı Eₒₒ için boyutsuz komşuluk korelasyonunu aşağıdaki gibi tanımlıyoruz.

$$
C_1 = \frac{E_{\mathrm{OO}}/E-f^2}{f(1-f)}.
$$

Pozitif değerler OH ve F'nin aynı türden siteyle komşuluğunun artmasına karşılık gelir. Bu ölçü tek komşuluk kabuğunu özetler. Aynı C₁ ile geometrik boyutları ve biçimleri farklı kümeler elde edilebilir. Kendine dönen bağ bulunmayan sonlu bir yüzde, sabit bileşimle rastgele atama ortalaması −1/(Nₛ−1) olur. `fixed_composition_random_mean` küçük hücrelerde aynı site görüntülerini de hesaba katar. Yalnızca OH veya yalnızca F içeren bir yüzün korelasyonu tanımsızdır.

`pattern: "correlated"` seçilince `surface_correlation_top` ve `surface_correlation_bottom` ile her yüzün hedefini girin. Her iki tabakanın yüzleri ayrı ayrı örneklenir. Bütün olası site atamaları en fazla 24.000 ise üretici gerekirse hepsini değerlendirir. Büyük hücrede en çok 40.000 sabit bileşimli OH↔F takası dener. `face_occupancy[tabaka][yüz].neighbor_statistics` hedefi, gerçekleşen değeri, farkı ve sayılan bağları verir. `exact` hedefin gerçekleştiğini, `nearest_feasible_proven` en yakın değerin ispatlandığını, `best_found_bounded_search` tarama bütçesinde bulunan en iyi değeri gösterir. Simülasyon koşulları karşılaştırılırken gerçekleşen değerleri kullanın. OH/F haritaları geometrik model girdileridir. Gerçek Ti₃C₂Tₓ için yüzey sonlanmalarının yerel çevresini inceleyen katı hâl NMR kaynağı [Hope ve arkadaşları](https://doi.org/10.1039/C6CP00330C).

İsteğe bağlı `vacuum_padding_A` üstte 0–100 Å boşluk bırakır. Varsayılan 0 Å eski hücre boyutunu korur. Bu boşluk başlangıçtaki çözücü hacmi, su/iyon sayısı ve koordinatlardan bağımsızdır. Üst z sınırı açık pilotlarda parçacık kaybını ayrıca izleyin.

## Yeniden üretim

`sequence_seed`, `tacticity_seed`, `conformer_seed` ve `surface_seed` ayrı değiştirilebilir. Yeni kontrollerden birini kullandığınızda verilmemiş tohumlar ana tohum ve bileşen adının SHA-256 özetinden türetilir. Aynı ana tohumlu eski reçeteler eski rastgelelik akışlarını korur. `seed_streams` gerçek değerleri saklar. Yalnızca dizi veya yüzey düzeni etkisini sınarken eşlenmiş koşulları kullanın. Bağımsız yapı tekrarları için ana tohumu değiştirin. Paketleme aşamasında yerleştirme tohumu ayrıca seçilebilir.

[Örnek reçete](../examples/studio/controlled_interface.json) ile iki tabakalı hidratasyon geometrisi oluşturulabilir. Altı yapının üretimini yeniden çalıştıran komut aşağıdadır.

```bash
python scripts/verify_studio_microstructure.py --output artifacts/interface_microstructures.json
```

[Ölçüm kaydı](../evidence/studio/interface_microstructures.json) atom sayısı, formül, iki üst yüzde gerçekleşen korelasyon, kesin asetat konumları ve koordinatların SHA-256 özetlerini içerir. Bu hesap geometrik ve kimyasal yapı denetimidir. Adsorpsiyon ve dispersiyon ölçümleri ayrıca doğrulama gerektirir. Atomların yükleri, kuvvet alanı ve dinamik kapsamı topoloji adımının çıktısında belgelenir.
