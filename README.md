# ⚙️ DualSPHysics_HIP
**Versión adaptada de DualSPHysics v5.0 al entorno HIP/ROCm**

Este repositorio corresponde a la **traducción y adaptación progresiva** del código base de DualSPHysics 5.0 optimizado en CUDA, para su ejecución sobre **GPUs AMD**, utilizando la plataforma **ROCm (Radeon Open Compute)** y el framework **HIP (Heterogeneous-compute Interface for Portability)**.  
Proyecto disponible en: [github.com/Nicolasg2c/DualSPHysics](https://github.com/Nicolasg2c/DualSPHysics)

---

## 🧩 Estructura del repositorio

```
DualSPHysics_HIP/
├── bin/                    # Binarios compilados
├── doc/                    # Documentación general y técnica
├── examples/               # Casos de prueba y ejemplos de simulación
├── src/                    # Código fuente principal adaptado a HIP
│   ├── src_CUDA/           # src adaptado a HIP
│   ├── src_HIPv01/         # src 1ra versión adaptada a HIP
|   ├── src_HIPv02/         # src 2da versión adaptada a HIP
│   └── src_HIPv03/         # src 3ra versión adaptada a HIP
├─ src01/VS
├─ src_extra/ToVTK_v5/
├─ src_mphase/
├── CHANGES.txt             # Registro de cambios y versiones
├── Files_DualSPHysics_v5.0.pdf  # Documentación original de referencia
├── LICENSE                 # Licencia LGPL-2.1
├── chpermissions.sh        # Script auxiliar para permisos de ejecución
└── README.md               
```

---

## 🧾 Licencia

Este proyecto se distribuye bajo la **GNU Lesser General Public License v2.1 (LGPL-2.1)**.  
Consulta el archivo `LICENSE` para más detalles.

---

## 👥 Autores

Proyecto desarrollado por:  
**Wilmer Farfán** y **Fabián Sánchez**  

Como parte del trabajo de grado titulado:  
> *“Análisis de la portabilidad de la implementación de métodos numéricos de hidrodinámica de partículas suaves en diferentes plataformas y frameworks CPU/GPU.”*

