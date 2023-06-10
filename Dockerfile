FROM mcr.microsoft.com/dotnet/runtime:6.0 AS base
WORKDIR /app

FROM mcr.microsoft.com/dotnet/sdk:6.0 AS build
WORKDIR /src
COPY ["Lab 1.csproj", "./"]
RUN dotnet restore "Lab 1.csproj"
COPY . .
WORKDIR "/src/"
RUN dotnet build "Lab 1.csproj" -c Release -o /app/build

FROM build AS publish
RUN dotnet publish "Lab 1.csproj" -c Release -o /app/publish

FROM base AS final
WORKDIR /app
COPY --from=publish /app/publish .
ENTRYPOINT ["dotnet", "Lab 1.dll"]
