%% Inicialización
clear all; close all; clc;
spm('defaults', 'pet');
spm_jobman('initcfg');

%% Función auxiliar para limpiar temp
function cleanTempDir(temp_dir)
   % Esperar un momento
   pause(2);
   
   % Intentar eliminar archivos primero
   if exist(temp_dir, 'dir')
       files = dir(fullfile(temp_dir, '*'));
       for k = 1:length(files)
           if ~strcmp(files(k).name,'.') && ~strcmp(files(k).name,'..')
               try
                   delete(fullfile(temp_dir, files(k).name));
               catch
                   % Si no se puede eliminar, ignorar y continuar
               end
           end
       end
       
       % Ahora intentar eliminar el directorio
       try
           rmdir(temp_dir, 's');
           mkdir(temp_dir);
       catch
           % Si no se puede eliminar, al menos está vacío
       end
   end
end

%% Definir directorios
% Directorio base donde están las imágenes
base_dir = 'C:/Users/juana/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM';
img_dir = fullfile(base_dir, 'PETimg_masked for simulations');
spm_dir = fullfile(base_dir, 'Results/z35/1vsGroup/SPM');

% Verificar directorios
if ~exist(img_dir, 'dir')
   error('Directorio de imágenes no encontrado');
end
if ~exist(spm_dir, 'dir')
   mkdir(spm_dir);
   disp('Directorio SPM creado');
end

%% Identificar archivos control (w00)
cd(img_dir);
control_files = dir('*w00*.nii');
num_controls = length(control_files);

% Verificar que encontramos los controles
if num_controls == 0
   error('No se encontraron archivos de control (w00)');
else
   disp(['Se encontraron ' num2str(num_controls) ' archivos de control']);
end

%% Identificar archivos de pacientes (excluyendo w79, w413 y niveles 2,6)
regiones = {'roiAD', 'w214', 'w271', 'w32'};
niveles = {'1', '4', '8'};
patient_files = [];

% Construir lista de archivos válidos
for r = 1:length(regiones)
   for n = 1:length(niveles)
       pattern = ['*_' regiones{r} '_0_' niveles{n} '_*.nii'];
       files = dir(pattern);
       patient_files = [patient_files; files];
   end
end

num_patients = length(patient_files);
if num_patients == 0
   error('No se encontraron archivos de pacientes');
else
   disp(['Se encontraron ' num2str(num_patients) ' archivos de pacientes']);
end

%% Crear carpeta temporal
temp_dir = fullfile(base_dir, 'Results/z35/1vsGroup/SPM/temp');
if ~exist(temp_dir, 'dir')
   mkdir(temp_dir);
end

%% Preparar rutas de controles (fuera del bucle principal)
control_paths = cell(num_controls, 1);
for i = 1:num_controls
   control_paths{i} = fullfile(img_dir, control_files(i).name);
end

%% Procesamiento de todos los pacientes
% Archivo de log para errores
log_file = fullfile(base_dir, 'Results/z35/1vsGroup/SPM/proceso_log.txt');
fid = fopen(log_file, 'a');
fprintf(fid, '\nNueva ejecución: %s\n', datestr(now));

% Bucle principal
for i = 1:length(patient_files)
   current_patient = patient_files(i).name;
   
   % Extraer información del paciente
   parts = split(current_patient, '_');
   control_num = parts{2};
   region = parts{4};
   level = parts{6};
   new_name = sprintf('binary_%s_%s_%s.nii', control_num, region, level);
   dest_file = fullfile(spm_dir, new_name);
   
   % Verificar si ya existe
   if exist(dest_file, 'file')
       msg = sprintf('Saltando paciente %d/%d: %s (ya existe)\n', ...
                    i, length(patient_files), new_name);
       fprintf(msg);
       fprintf(fid, msg);
       continue;
   end
   
   try
       msg = sprintf('Procesando paciente %d/%d: %s\n', ...
                    i, length(patient_files), current_patient);
       fprintf(msg);
       fprintf(fid, msg);
       
       % Limpiar matlabbatch al inicio de cada iteración
       clear matlabbatch
       
       % Preparar rutas
       test_patient = fullfile(img_dir, current_patient);
       
       % Configurar matlabbatch para SPM
       matlabbatch{1}.spm.stats.factorial_design.dir = {temp_dir};
       matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = {test_patient};
       matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = control_paths;
       matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
       matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
       matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
       matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
       matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
       matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
       matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
       matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
       matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
       
       % Ejecutar el diseño factorial
       spm_jobman('run', matlabbatch);
       
       % Configurar la estimación
       clear matlabbatch
       matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(temp_dir, 'SPM.mat')};
       matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
       matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
       spm_jobman('run', matlabbatch);
       
       % Configurar el contraste [-1 1]
       clear matlabbatch
       matlabbatch{1}.spm.stats.con.spmmat = {fullfile(temp_dir, 'SPM.mat')};
       matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Patient < Control';
       matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [-1 1];
       matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
       spm_jobman('run', matlabbatch);
       
       % Configurar thresholding
       clear matlabbatch
       matlabbatch{1}.spm.stats.results.spmmat = {fullfile(temp_dir, 'SPM.mat')};
       matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
       matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
       matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
       matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
       matlabbatch{1}.spm.stats.results.conspec.extent = 0;
       matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
       matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
       spm_jobman('run', matlabbatch);

       % Convertir a binary
       spmT = spm_select('FPList', temp_dir, '^spmT.*\.nii$');
       binary_file = fullfile(temp_dir, 'binary.nii');
       if exist(spmT, 'file')
           V = spm_vol(spmT);
           Y = spm_read_vols(V);
           Y = Y > 0;  % Convertir a binario
           V.fname = binary_file;
           spm_write_vol(V, Y);
           
           % Mover y renombrar
           movefile(binary_file, dest_file);
           
           % Limpiar carpeta temporal
           cleanTempDir(temp_dir);
           
           msg = sprintf('Completado con éxito: %s\n', new_name);
           fprintf(msg);
           fprintf(fid, msg);
       else
           error('No se encontró el archivo spmT');
       end
       
       % Pausa entre iteraciones
       pause(1);
       
   catch ME
       % Registrar error
       msg = sprintf('ERROR en paciente %s: %s\n', ...
                    strrep(current_patient, '\', '\\'), strrep(ME.message, '\', '\\'));
       fprintf(msg);
       fprintf(fid, msg);
       
       % Limpiar temp si hubo error
       cleanTempDir(temp_dir);
   end
end

% Cerrar log
fclose(fid);