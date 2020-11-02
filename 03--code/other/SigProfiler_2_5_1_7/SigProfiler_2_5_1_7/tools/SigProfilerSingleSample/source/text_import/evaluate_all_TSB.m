function output = evaluate_all_TSB(input)

    totalSamples = size(input.originalGenomes,2);

    %% Pre-generating output matrices
    output.C_to_A_p = zeros(totalSamples, 1);      output.C_to_A_d = zeros(totalSamples, 1);
    output.C_to_G_p = zeros(totalSamples, 1);      output.C_to_G_d = zeros(totalSamples, 1);
    output.C_to_T_p = zeros(totalSamples, 1);      output.C_to_T_d = zeros(totalSamples, 1);
    output.T_to_A_p = zeros(totalSamples, 1);      output.T_to_A_d = zeros(totalSamples, 1);
    output.T_to_C_p = zeros(totalSamples, 1);      output.T_to_C_d = zeros(totalSamples, 1);
    output.T_to_C_ATN_p = zeros(totalSamples, 1);  output.T_to_C_ATN_d = zeros(totalSamples, 1);
    output.T_to_G_p = zeros(totalSamples, 1);      output.T_to_G_d = zeros(totalSamples, 1); 

    % C > A
    C_to_A_U = strcmp('C>A',input.types) & strcmp('U', input.strandTypes);
    C_to_A_T = strcmp('C>A',input.types) & strcmp('T', input.strandTypes);

    % C > G
    C_to_G_U = strcmp('C>G',input.types) & strcmp('U', input.strandTypes);
    C_to_G_T = strcmp('C>G',input.types) & strcmp('T', input.strandTypes);

    % C > T
    C_to_T_U = strcmp('C>T',input.types) & strcmp('U', input.strandTypes);
    C_to_T_T = strcmp('C>T',input.types) & strcmp('T', input.strandTypes);

    % T > A
    T_to_A_U = strcmp('T>A',input.types) & strcmp('U', input.strandTypes);
    T_to_A_T = strcmp('T>A',input.types) & strcmp('T', input.strandTypes);

    % T > C
    T_to_C_U = strcmp('T>C',input.types) & strcmp('U', input.strandTypes);
    T_to_C_T = strcmp('T>C',input.types) & strcmp('T', input.strandTypes);

    % T > C @ ApTpN
    T_to_C_U_ATN = strcmp('T>C',input.types) & strcmp('U', input.strandTypes) & ...
                   ( strcmp('ATA', input.subtypes) | strcmp('ATC', input.subtypes) ...
                   | strcmp('ATG', input.subtypes)  | strcmp('ATT', input.subtypes) );
    T_to_C_T_ATN = strcmp('T>C',input.types) & strcmp('T', input.strandTypes) & ...
                   ( strcmp('ATA', input.subtypes) | strcmp('ATC', input.subtypes) ...
                   | strcmp('ATG', input.subtypes)  | strcmp('ATT', input.subtypes) );

    % T > G
    T_to_G_U = strcmp('T>G',input.types) & strcmp('U', input.strandTypes);
    T_to_G_T = strcmp('T>G',input.types) & strcmp('T', input.strandTypes);

    for i = 1 : totalSamples
        genome = input.originalGenomes(:,i);

        % Testing strand bias using Fisher exact test
        [output.C_to_A_p(i), output.C_to_A_d(i)] = evaluateTSB(C_to_A_U, C_to_A_T, genome);
        [output.C_to_G_p(i), output.C_to_G_d(i)] = evaluateTSB(C_to_G_U, C_to_G_T, genome);
        [output.C_to_T_p(i), output.C_to_T_d(i)] = evaluateTSB(C_to_T_U, C_to_T_T, genome);
        [output.T_to_A_p(i), output.T_to_A_d(i)] = evaluateTSB(T_to_A_U, T_to_A_T, genome);
        [output.T_to_C_p(i), output.T_to_C_d(i)] = evaluateTSB(T_to_C_U, T_to_C_T, genome);
        [output.T_to_C_ATN_p(i), output.T_to_C_ATN_d(i)] = evaluateTSB(T_to_C_U_ATN, T_to_C_T_ATN, genome);
        [output.T_to_G_p(i), output.T_to_G_d(i)] = evaluateTSB(T_to_G_U, T_to_G_T, genome);
    end
    
end

