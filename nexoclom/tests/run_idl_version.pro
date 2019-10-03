pro run_idl_version, inputfile, overwrite=overwrite

if n_elements(overwrite) EQ 0 then overwrite = 0

inputs = inputs_restore(inputfile)
sp = inputs.options.atom
if sp EQ 'Ca' $
    then orbit = 1489 $
    else orbit = 1541

data = load_mascs_data(sp, orbit)
inputs.geometry.taa = median(*data.taa)

modeldriver, inputs, 1e5, overwrite=overwrite
outputfile = modeloutput_search(inputs)

formatfile0 = 'inputfiles/MercuryEmissionIDL'+sp+'.format'
image = produce_results(inputs, formatfile0)
im = *image.image

formatfile1 = 'inputfiles/MercuryColumnIDL.format'
column = produce_results(inputs, formatfile1)
col = *column.image

params = list('MESSENGER', 'intensity', orbit, 3*!dtor, sp)
format_inten = make_format_structure(params)
model = produce_results(inputs, format_inten, data=data)
radiance = *model.radiance
packets = *model.npackets

savefile = 'idloutputs/' + file_basename(inputfile) + '.sav'

save, file=savefile, outputfile, im, col, radiance, packets

print, 'Completed Successfully'

end
